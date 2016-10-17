# -*- coding: utf-8 -*-
"""
Helper functions to plot results.
"""
import nilearn.image      as niimg
import numpy              as np
import os.path as op
import pandas             as pd
import re
import scipy.io           as sio
from   boyle.nifti.utils  import filter_icc
from   nilearn.image      import iter_img
from   nilearn.masking    import apply_mask

from .ica_loadings import (filter_ics,
                           get_largest_blobs,
                           build_raw_loadings_table,
                           add_groups_to_loadings_table,
                           )
from ..utils.files import fetch_one_file
from ..interfaces.nilearn.plot import (plot_all_components,
                                       plot_ica_components,
                                       plot_multi_slices,
                                       plot_overlays)


def plot_connectivity_matrix(connectivity_matrix, label_names):
    """ Plot the connectivity matrix. """
    # Display the correlation matrix
    import numpy as np
    from matplotlib import pyplot as plt

    plt.figure(figsize=(10, 10))

    # Mask out the major diagonal
    np.fill_diagonal(connectivity_matrix, 0)
    fig = plt.imshow(connectivity_matrix, interpolation="nearest", cmap="RdBu_r",
                     vmax=0.8, vmin=-0.8)
    plt.colorbar()
    # And display the labels
    _ = plt.xticks(range(len(label_names)), label_names, rotation=90)
    _ = plt.yticks(range(len(label_names)), label_names)

    return fig


def plot_ica_results(ica_result_dir, application, mask_file='', zscore=2., mode='+', bg_img=''):
    """  Use nilearn to plot results from different ICA analysis tools, given the ICA result folder path.

    Parameters
    ----------
    ica_result_dir: str
        Path to the ICA output folder or the ICA components volume file.

    application: str
        Choicese: ('nilearn', 'sbm', 'gift')

    mask_file: str
        Path to the brain mask file to be used for thresholding.

    mode: str
        Choices: '+' for positive threshold,
                 '+-' for positive and negative threshold and
                 '-' for negative threshold.

    zscore: int
        Value of the Z-score thresholding.

    bg_img: str
        Path to a background image.
        If empty will use the SPM canonical brain image at 2mm.
    """
    if application == 'sbm': # SBM ICA (3D sources)
        plotter = SBMICAResultsPlotter(ica_result_dir)
    elif application == 'gift': # group GIFT ICA (4D sources)
        plotter = GIFTICAResultsPlotter(ica_result_dir)
    elif application == 'canica': # group ICA (4D sources)
        plotter = CanICAResultsPlotter(ica_result_dir)
    else:
        raise NotImplementedError('Got no ICAResultsPlotter for application {}.'.format(application))

    plotter.fit(mask_file=mask_file, mode=mode, zscore=zscore)
    plotter.plot_icmaps(bg_img=bg_img)


class ICAResultsPlotter(object):
    """ Use nilearn to plot results from different ICA analysis tools, given the ICA result folder path.

    Parameters
    ----------
    ica_result_dir: str
        Path to the ICA output folder or the ICA components volume file.

    mask_file: str
        Path to the brain mask file to be used for thresholding.

    mode: str
        Choices: '+' for positive threshold,
                 '+-' for positive and negative threshold and
                 '-' for negative threshold.

    zscore: int or float
        Value of the Z-score thresholding.
    """
    _comps_fname = None

    def __init__(self, ica_result_dir):
        if not op.exists(ica_result_dir):
            raise IOError('Expected an existing file or folder, but could not find {}.'.format(ica_result_dir))

        self.ica_dir = op.expanduser(ica_result_dir)

    def fit(self,  mask_file='', mode='+', zscore=2):
        """ Process/filter/threshold the output to make it ready for plot.

        Parameters
        ----------
        mask_file: str
            Path to the brain mask file to be used for thresholding.

        mode: str
            Choices: '+' for positive threshold,
                     '+-' for positive and negative threshold and
                     '-' for negative threshold.

        zscore: int or float
            Value of the Z-score thresholding.
        """
        self.mask_file = mask_file
        self.mode = mode
        self.zscore = zscore
        self._icc_imgs = None
        self._update(force=True)

    def _fetch_components_file(self):
        if self._comps_fname is None:
            raise NotImplementedError('This is a generic class to support the output from different applications,'
                                      'please use a derived class from ICAResultsPlotter.')

        return fetch_one_file(self.ica_dir, self._comps_fname)

    def _filter_ic_imgs(self, ic_file):
        # filter the ICC if mask and threshold are set
        if self.mask_file and self.zscore > 0:
            mask = niimg.load_img(self.mask_file)
            icc_imgs = [filter_icc(icc, mask=mask, thr=self.zscore, zscore=True, mode=self.mode)
                        for icc in iter_img(ic_file)]
        else:
            icc_imgs = [niimg.load_img(ic) for ic in ic_file]
        return icc_imgs

    def _check_output(self):
        raise NotImplementedError('This is a generic class to support the output from different applications,'
                                  'please use a derived class from ICAResultsPlotter.')

    def _update(self, force=False):
        if not hasattr(self, '_icc_imgs'):
            raise RuntimeError('You should call `fit()` before using this function.')

        if self._icc_imgs is not None and not force:
            return

        self._check_output()

        ic_file = self._fetch_components_file()
        self._icc_imgs = list(filter_ics(ic_file, mask=self.mask_file, zscore=self.zscore))

    def plot_icmaps(self, outtype='png', **kwargs):
        """ Plot the thresholded IC spatial maps and store the outputs in the ICA results folder.
        Parameters
        ----------
        outtype: str
            Extension (without the '.') of the output files, will specify which plot image file you want.

        Returns
        -------
        all_icc_plot_f: str

        iccs_plot_f: str

        sliced_ic_plots: list of str
        """
        # specify the file paths
        all_icc_plot_f  = op.join(self.ica_dir, 'all_components_zscore_{}.{}'.format(self.zscore, outtype))
        iccs_plot_f     = op.join(self.ica_dir,  'ic_components_zscore_{}.{}'.format(self.zscore, outtype))
        icc_multi_slice = op.join(self.ica_dir, 'ic_map_{}_zscore_{}.{}')

        # make the plots
        fig1 = plot_ica_components(self._icc_imgs, **kwargs)
        fig1.savefig(iccs_plot_f, facecolor=fig1.get_facecolor(), edgecolor='none')

        fig2 = plot_all_components(self._icc_imgs, **kwargs)
        fig2.savefig(all_icc_plot_f, facecolor=fig2.get_facecolor(), edgecolor='none')

        # make the multi sliced IC plots
        sliced_ic_plots = []
        for i, img in enumerate(iter_img(self._icc_imgs)):
            fig3 = plot_multi_slices(img,
                                     cut_dir="z",
                                     n_cuts=24,
                                     n_cols=4,
                                     title="IC {}\n(z-score {})".format(i+1, self.zscore),
                                     title_fontsize=32,
                                     plot_func=None,
                                     **kwargs)

            # prepare the output file name/path
            out_f = icc_multi_slice.format(i+1, self.zscore, outtype)
            fig3.savefig(out_f, facecolor=fig3.get_facecolor(), edgecolor='none')
            sliced_ic_plots.append(out_f)

        return all_icc_plot_f, iccs_plot_f, sliced_ic_plots


class CanICAResultsPlotter(ICAResultsPlotter):
    """ Use nilearn to plot results from CanICA and DictLearing tools from nilearn, given the ICA result folder path.

    Parameters
    ----------
    ica_result_dir: str
        Path to the ICA output folder or the ICA components volume file.
    """
    _comps_fname = 'canica_resting_state.nii.gz'


class MIALABICAResultsPlotter(ICAResultsPlotter):
    """ Use nilearn to plot results from MIALAB (GIFT and SBM) tools, given the ICA result folder path.

    Parameters
    ----------
    ica_result_dir: str
        Path to the ICA output folder or the ICA components volume file.
    """
    # define the input file patterns
    _comps_fname     = '*_component_ica*'
    _loadings_fname  = '*_loading_coeff*'
    _subjects_fname  = '*Subject.mat'
    _mask_fname      = '*Mask.hdr'

    def __init__(self, ica_result_dir):
        super(MIALABICAResultsPlotter, self).__init__(ica_result_dir=ica_result_dir)
        self._subjid_pat = ''

    def _parse_groups_file(self, group_labels_file):
        """ Return the list of groups and subject ids from the group_labels_file and the Subjects.mat file.
        the file paths inside the *Subject.mat file.
        One file is `subject_loadings.xls` which has the loadings as is, with the subjects IDs and group.
        The other file is `subject_group_loadings.xls` which has the loading signs changed according to
        the average correlation value of the "main" region of each of the IC spatial maps.

        Parameters
        ----------
        group_labels_file: str
            A CSV file with two columns: "subject_id" and "group".
            The subject_ids must be in the paths contained in the Subject.mat
            file and match the `subjid_pat` argument.

        Returns
        -------
        groups: pandas.DataFrame
        """
        groups = None
        if group_labels_file and op.exists(group_labels_file):
            groups = pd.read_csv(group_labels_file)
            if not 'subject_id' in groups.columns or \
               not 'group' in groups.columns:
                raise AttributeError("Please add columns names 'subject_id' "
                                     "and 'group' to the group labels file.")
        return groups

    def _load_components(self):
        """ Return components image file and check shape match."""
        compsf = self._fetch_components_file()
        comps_img = niimg.load_img(compsf)
        return comps_img

    def _fetch_components_file(self):
        return fetch_one_file(self.ica_dir, self._comps_fname,
                              file_extension='.nii',
                              extra_prefix='*mean',
                              extra_suffix='_s_all*')

    def _get_subject_files(self):
        """ Load the .mat file with subjects lists, mainly to get the order in
        which subjects are introduced in the other matrices.
        """
        from itertools import chain

        subjsf   = fetch_one_file(self.ica_dir, self._subjects_fname)
        mat_file = sio.loadmat(subjsf)['files']
        return [f.strip() for f in list(chain.from_iterable(chain.from_iterable(chain.from_iterable(mat_file))))]

    def _get_subject_ids(self, subjid_pat):
        """ Return the list of subject ids parsed from the file paths present in the Subjects.mat file.

        Parameters
        ----------
        subjid_pat: regext str
            A search regex pattern that returns one group element that
            contains the subject id.
            This will be used to search for subject_id in the file paths
            contained in the Subjects.mat file.

        Returns
        -------
        patids: list[str]
        """
        self._subjid_pat = subjid_pat

        subj_files = self._get_subject_files()

        # process the paths to extract the patient IDs given by `patid_re`
        in_files = [f.split(',')[0] for f in subj_files]
        patid_re = re.compile(subjid_pat, re.I)
        patids   = [patid_re.search(f).group() for f in in_files]

        return patids

    def fit(self,  mask_file='', mode='+', zscore=2):
        """ Process/filter/threshold the output to make it ready for plot.
        If no mask_file is set, will pick the output from GIFT.

        Parameters
        ----------
        mask_file: str
            Path to the brain mask file to be used for thresholding.

        mode: str
            Choices: '+' for positive threshold,
                     '+-' for positive and negative threshold and
                     '-' for negative threshold.

        zscore: int or float
            Value of the Z-score thresholding.
        """
        if not mask_file:
            mask_file = fetch_one_file(self.ica_dir, self._mask_fname)

        super(MIALABICAResultsPlotter, self).fit(mask_file=mask_file,
                                                 mode=mode,
                                                 zscore=zscore)


class GIFTICAResultsPlotter(MIALABICAResultsPlotter):
    """ Use nilearn to plot results from the MIALAB SBM tool, given the ICA result folder path.

    Parameters
    ----------
    ica_result_dir: str
        Path to the ICA output folder or the ICA components volume file.
    """
    _tcs_fname = '*_timecourses_ica*.nii'

    def _check_output(self):
        comps_img = self._load_components()
        tcs_img   = self._load_timecourses()

        # read the groups file
        subj_files   = self._get_subject_files()
        n_timepoints = len(subj_files)

        n_ics = tcs_img.shape[1]
        if comps_img.get_data().shape[-1] != n_ics:
            raise AttributeError('Shape mismatch between timecourses matrix {} '
                                 'and components data shape {}.'.format(tcs_img.shape,
                                                                        comps_img.get_data().shape))

        # check shapes
        #if n_timepoints != tcs_img.shape[0]:
        #    raise AttributeError('Shape mismatch between list of subjects {} '
        #                         'and timecourses data shape {}.'.format(n_timepoints,
        #                                                                 tcs_img.shape))

    def _load_timecourses(self):
        """ Return the timecourses image file and checks if the shape is correct."""
        # load the timecourses file
        tcsf = fetch_one_file(self.ica_dir, self._tcs_fname, extra_prefix='*tmap')
        tcs = niimg.load_img(tcsf).get_data()
        return tcs


    def _calculate_goodness_of_fit(self):
        """ Return the goodness-of-fit values from a GIFT result."""
        pass
        #TODO


    def goodness_of_fit_df(self, group_labels_file, subjid_pat=r'(?P<patid>[a-z]{2}_[0-9]{6})'):
        """ Return a pandas.DataFrame ready for an excel file with:
        - the subject IDs taken from the file paths inside the *Subject.mat file,
        - the groups taken from `group_labels_file`, and
        - the goodness-of-fit measures for each subject and independent component.

        Parameters
        ----------
        group_labels_file: str
            A CSV file with two columns: "subject_id" and "group".
            The subject_ids must be in the paths contained in the Subject.mat
            file and match the `subjid_pat` argument.

        subjid_pat: regext str
            A search regex pattern that returns one group element that
            contains the subject id.
            This will be used to *search* for subject_id in the file paths
            contained in the Subjects.mat file.

        Returns
        -------
        gof_df: pandas.DataFrame
        """
        # make sure file exists
        if not op.exists(group_labels_file):
            raise FileNotFoundError('The file {} has not been found.'.format(group_labels_file))

        # make sure this object has been .fit()
        self._update()

        # read the groups file
        groups = self._parse_groups_file(group_labels_file=group_labels_file)
        patids = self._get_subject_ids(subjid_pat=subjid_pat)

        # calculate the goodness of fit
        gofs = self._calculate_goodness_of_fit(patids)

        # build the goodness-of-fit table
        df = build_raw_loadings_table(gofs, patids)
        df = add_groups_to_loadings_table(df, groups)

        return df


class SBMICAResultsPlotter(MIALABICAResultsPlotter):
    """ Use nilearn to plot results from the MIALAB SBM tool, given the ICA result folder path.

    Parameters
    ----------
    ica_result_dir: str
        Path to the ICA output folder or the ICA components volume file.
    """
    def _load_loadings(self):
        loadf = fetch_one_file(self.ica_dir, self._loadings_fname, file_extension='.nii')
        loads = niimg.load_img(loadf).get_data()
        return loads

    def _check_output(self):
        # read the groups file
        subj_files   = self._get_subject_files()
        n_timepoints = len(subj_files)

        # load the loadings file
        loads = self._load_loadings()

        # check shapes
        if n_timepoints != loads.shape[0]:
            raise AttributeError('Shape mismatch between list of subjects {} '
                                 'and loadings data shape {}.'.format(n_timepoints,
                                                                      loads.shape))

        # load components image file to check shape match
        # this doesn't make much sense to be here (because it does not use the components image
        # but it checks if the ICA output is consistent).
        compsf = fetch_one_file(self.ica_dir, self._comps_fname)
        comps_img = niimg.load_img(compsf)

        n_ics = loads.shape[1]
        if comps_img.get_data().shape[-1] != n_ics:
            raise AttributeError('Shape mismatch between loadings matrix {} '
                                 'and components data shape {}.'.format(loads.shape,
                                                                        comps_img.get_data().shape))

    def simple_loadings_df(self, group_labels_file, subjid_pat=r'(?P<patid>[a-z]{2}_[0-9]{6})'):
        """ Return a pandas.DataFrame spreadsheet ready for an excel file with the subject IDs taken from
        the file paths inside the *Subject.mat file.
        One file is `subject_loadings.xls` which has the loadings as is, with the subjects IDs and group.
        The other file is `subject_group_loadings.xls` which has the loading signs changed according to
        the average correlation value of the "main" region of each of the IC spatial maps.

        Parameters
        ----------
        group_labels_file: str
            A CSV file with two columns: "subject_id" and "group".
            The subject_ids must be in the paths contained in the Subject.mat
            file and match the `subjid_pat` argument.

        subjid_pat: regext str
            A search regex pattern that returns one group element that
            contains the subject id.
            This will be used to search for subject_id in the file paths
            contained in the Subjects.mat file.

        Returns
        -------
        loadings_df: pandas.DataFrame
        """
        # make sure file exists
        if not op.exists(group_labels_file):
            raise FileNotFoundError('The file {} has not been found.'.format(group_labels_file))

        # make sure this object has been .fit()
        self._update()

        # read the groups file
        groups = self._parse_groups_file(group_labels_file=group_labels_file)
        patids = self._get_subject_ids(subjid_pat=subjid_pat)

        # load the loadings file
        loads = self._load_loadings()

        # build the raw loadings table
        df = build_raw_loadings_table(loads, patids)
        df = add_groups_to_loadings_table(df, groups)
        #df.to_excel(op.join(ica_out_dir, rawloadings_filename))

        return df

    def weighted_loadings_df(self, group_labels_file, subjid_pat=r'(?P<patid>[a-z]{2}_[0-9]{6})'):
        """ Return a pandas.DataFrame with the loadings, subject IDs and groups.
        This version of the sheet has the loadings signs changed for each group of subjects.
        Multiply the largest blob of each IC spatial map (see `plot_ic_largest_blobs`) and the sign of the
        average loading of each group.

        Parameters
        ----------
        group_labels_file: str
            A CSV file with two columns: "subject_id" and "group".
            The subject_ids must be in the paths contained in the Subject.mat
            file and match the `subjid_pat` argument.

        subjid_pat: regext str
            A search regex pattern that returns one group element that
            contains the subject id.
            This will be used to search for subject_id in the file paths
            contained in the Subjects.mat file.

        Returns
        -------
        loadings_df: pandas.DataFrame
        """
        # make sure file exists
        if not op.exists(group_labels_file):
            raise FileNotFoundError('The file {} has not been found.'.format(group_labels_file))

        # make sure this object has been .fit()
        self._update()

        # let's first pick the simple version of the loadings
        df = self.simple_loadings_df(group_labels_file, subjid_pat=subjid_pat)

        # get the values of the largest blobs in the filtered IC maps
        blobs = get_largest_blobs(self._icc_imgs)

        # apply the blob mask to each corresponding IC map
        masks = [apply_mask(ic_map, blob) for ic_map, blob in zip(self._icc_imgs, blobs)]

        # calculate the avg per blob
        blob_avgs = [mask.mean() for mask in masks]

        # multiply the avg blob value with the group loadings
        blob_signs = np.sign(blob_avgs)
        n_ics      = len(blob_avgs)

        df[list(range(1, n_ics+1))] = df[list(range(1, n_ics+1))] * blob_signs

        # put the dataframe back into one whole ungrouped dataframe
        return df

    def plot_icmaps_and_blobs(self, outfile, bg_img=None, **kwargs):
        """ Plot the IC maps with the largest blobs overlaid.
        Parameters
        ----------
        outfile: str
            Path to output plot file.

        bg_img: str
            Path to image file.
        """
        # make sure this object has been .fit()
        self._update()
        fig = plot_overlays(self._icc_imgs, get_largest_blobs(self._icc_imgs),
                            bg_img=bg_img, figsize=(2.5, 3), **kwargs)

        # save fig
        fig.savefig(outfile, facecolor=fig.get_facecolor(), edgecolor='none')
