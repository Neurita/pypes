# -*- coding: utf-8 -*-
"""
Helper functions to plot results.
"""
import re
import os.path as op

import nilearn.image      as niimg
import numpy              as np
import pandas             as pd
import scipy.io           as sio
from   nilearn.input_data import NiftiMasker
from   nilearn.image      import iter_img
from   nilearn.masking    import apply_mask
from   nilearn._utils.niimg_conversions import check_niimg, _index_img
from   boyle.nifti.utils  import filter_icc

from .utils import (get_largest_blobs,
                    build_raw_loadings_table,
                    add_groups_to_loadings_table,)
from ..interfaces import (plot_all_components,
                          plot_ica_components,
                          plot_multi_slices,
                          plot_overlays)
from ..utils import fetch_one_file


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
        Choicese: ('nilearn', 'sbm', 'gift', 'gift-group')

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
    elif application == 'gift-group': # group GIFT ICA (4D sources)
        plotter = GIFTGroupICAResultsPlotter(ica_result_dir)
    elif application == 'gift': # single subject GIFT ICA (4D sources)
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

    def fit(self,  mask_file='', mode='+-', zscore=2):
        """ Process/filter/threshold the output to make it ready for plot.

        Parameters
        ----------
        mask_file: str
            Path to the brain mask file to be used for thresholding.

        mode: str
            Choices: '+' for positive threshold,
                     '+-' or '-+' for positive and negative threshold and
                     '-' for negative threshold.

        zscore: int or float
            Value of the Z-score thresholding.
        """
        self.mask_file = mask_file
        self.mode = mode
        self.zscore = zscore
        self._icc_imgs = None
        self._update(force=True)

    def _update(self, force=False):
        if self._icc_imgs is not None and not force:
            return

        self._check_output()
        self._icc_imgs = self._filter_ic_imgs(self._fetch_components_file())

    def is_fit(self):
        """ Return True if the self `fit` function has been called, False otherwise."""
        if not hasattr(self, '_icc_imgs'):
            return False
        else:
            return self._icc_imgs is not None

    def _check_is_fit(self):
        if not self.is_fit():
            raise RuntimeError('You should call `fit()` before using this function.')

    def _fetch_components_file(self):
        if self._comps_fname is None:
            raise NotImplementedError('This is a generic class to support the output from different applications,'
                                      'please use a derived class from ICAResultsPlotter.')

        return fetch_one_file(self.ica_dir, self._comps_fname)

    def _filter_ic_imgs(self, ic_file):
        if self.zscore > 0:
            do_zscore = True
        else:
            do_zscore = False

        mask = niimg.load_img(self.mask_file)
        return [filter_icc(icimg, mask=mask, thr=self.zscore, zscore=do_zscore, mode=self.mode)
                for icimg in iter_img(ic_file)]

    def _check_output(self):
        raise NotImplementedError('This is a generic class to support the output from different applications,'
                                  'please use a derived class from ICAResultsPlotter.')

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
    _tcs_fname       = '*_timecourses_ica*'
    _comps_fname     = '*_component_ica*'
    _subjects_fname  = '*Subject.mat'
    _mask_fname      = r'.*Mask.(hdr|nii)'

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

    def _load_loadings(self):
        loadf = fetch_one_file(self.ica_dir, self._tcs_fname)
        loads = niimg.load_img(loadf).get_data()
        return loads

    def _fetch_components_file(self):
        return fetch_one_file(self.ica_dir, self._comps_fname)

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

    @staticmethod
    def _index_img(img_file, index):
        """ Return the volume in `index` 4th dimension index of `img_file`.
        If the image is 3D and idx is zero will return the container of `img_file`.
        """
        imgs = check_niimg(img_file, ensure_ndim=4, atleast_4d=True)
        return _index_img(imgs, index)

    def load_mask(self):
        """ Return the mask image. """
        mask_file = fetch_one_file(self.ica_dir, self._mask_fname, pat_type='re.match')
        return niimg.load_img(mask_file)

    def load_subject_data(self, masked=False, **kwargs):
        """ Generator of the input data volumes in the order it was inserted in GIFT ICA.

        Returns
        -------
        subj_imgs: generator of niimg-like objects.
            The input data to the ICA.

        masked: bool
            If True will apply the mask to the data.

        kwargs: keyword arguments
            Keyword arguments for the NiftiMasker used to mask the data.
            The mask used is the one used in the `fit` function or specified
            in `_mask_fname`.

        Note
        ----
        This will actually read the paths in the Subjects.mat file. If you have moved
        or erased those files, this will not work as expected.
        """
        if masked:
            sessions = self._input_data_sessions()
            data = self._apply_mask_to_imgs(self._load_subject_data(), sessions=sessions, **kwargs)
        else:
            data = self._load_subject_data()

        return np.array(list(data))

    def load_components_data(self, masked=False, **kwargs):
        if masked:
            comps = self._load_components()
            n_comps = comps.shape[-1]
            return self._apply_mask_to_4dimg(comps, **kwargs)
        else:
            return self._load_components()

    def _input_data_sessions(self):
        """ Return a 1D array with a session number for each session
        in the input data.
        This is used for NiftiMasker as a value for the `sessions`
        parameter.

        Example
        -------
        Given input data:
        ["subject01.nii,1",
         "subject01.nii,2",
         "subject02.nii,1",]

        Would return: [0, 0, 1]
        """
        # get the subject file list
        subj_files = self._get_subject_files()

        # remove the numbers after the comma in this list
        files = np.array([f.split(',')[0] for f in subj_files])

        # make a dict that assign a number to each unique file
        file_int = {f: i for i, f in enumerate(np.unique(files))}

        # return the corresponding numbers for each file in the cleaned file list
        return np.array([file_int[f] for f in files])

    def _load_subject_data(self):
        img_file_list = self._get_subject_files()
        for line in img_file_list:
            img_file, idx = line.split(',')
            yield self._index_img(img_file, int(idx) - 1)

    def _apply_mask_to_img(self, img, **kwargs):
        masker = NiftiMasker(mask_img=self.load_mask(), **kwargs)
        return masker.fit_transform(img)

    def _apply_mask_to_4dimg(self, imgs, **kwargs):
        masker = NiftiMasker(mask_img=self.load_mask(), **kwargs)
        return (masker.fit_transform(img) for img in iter_img(imgs))

    def _apply_mask_to_imgs(self, imgs, **kwargs):
        masker = NiftiMasker(mask_img=self.load_mask(), **kwargs)
        return (masker.fit_transform(img) for img in imgs)

    def fit(self,  mask_file='', mode='+', zscore=2):
        """ Process/filter/threshold the output to make it ready for plot.
        If no mask_file is set, will pick the output from GIFT.

        Parameters
        ----------
        mask_file: str, optional
            Path to the brain mask file to be used for thresholding.

        mode: str, optional
            Choices: '+' for positive threshold,
                     '+-' for positive and negative threshold and
                     '-' for negative threshold.

        zscore: int or float, optional
            Value of the Z-score thresholding.
        """
        if not mask_file:
            mask_file = fetch_one_file(self.ica_dir, self._mask_fname, pat_type='re.match')

        super(MIALABICAResultsPlotter, self).fit(mask_file=mask_file,
                                                 mode=mode,
                                                 zscore=zscore)

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


class GIFTICAResultsPlotter(MIALABICAResultsPlotter):
    """ Use nilearn to plot results from the MIALAB ICA tool used to analyse one subject,
    given the ICA result folder path.

    Parameters
    ----------
    ica_result_dir: str
        Path to the ICA output folder or the ICA components volume file.
    """
    _tcs_fname   = r'.*_sub[0]+1_timecourses_ica_s1_.nii'
    _comps_fname = r'.*_sub[0]+1_component_ica_s1_.nii'

    def _check_output(self):
        comps_img = self._load_components()
        tcs_img   = self._load_timecourses()

        n_ics = tcs_img.shape[1]
        if comps_img.get_data().shape[-1] != n_ics:
            raise AttributeError('Shape mismatch between timecourses matrix {} '
                                 'and components data shape {}.'.format(tcs_img.shape,
                                                                        comps_img.get_data().shape))

        # read the groups file
        subj_files   = self._get_subject_files()
        n_timepoints = len(subj_files)

        # check shapes
        if n_timepoints != tcs_img.shape[0]:
            raise AttributeError('Shape mismatch between list of subjects {} '
                                 'and timecourses data shape {}.'.format(n_timepoints,
                                                                         tcs_img.shape))

    def _load_timecourses(self):
        """ Return the timecourses image file and checks if the shape is correct."""
        # load the timecourses file
        tcsf = fetch_one_file(self.ica_dir, self._tcs_fname, pat_type='re.match')
        tcs = niimg.load_img(tcsf).get_data()
        return tcs

    def _fetch_components_file(self):
        return fetch_one_file(self.ica_dir, self._comps_fname, pat_type='re.match')


class GIFTGroupICAResultsPlotter(MIALABICAResultsPlotter):
    """ Use nilearn to plot results from the MIALAB Group ICA tool used to analyse a group of subjects,
    given the ICA result folder path.

    Parameters
    ----------
    ica_result_dir: str
        Path to the ICA output folder or the ICA components volume file.
    """
    _tcs_fname   = '*_agg__timecourses_ica_.nii'
    _comps_fname = '*_agg__component_ica_.nii'

    def _check_output(self):
        comps_img = self._load_components()
        tcs_img   = self._load_timecourses()

        n_ics = tcs_img.shape[1]
        if comps_img.get_data().shape[-1] != n_ics:
            raise AttributeError('Shape mismatch between timecourses matrix {} '
                                 'and components data shape {}.'.format(tcs_img.shape,
                                                                        comps_img.get_data().shape))

        # read the groups file
        # subj_files   = self._get_subject_files()
        # n_timepoints = len(subj_files)

        # check shapes
        #if n_timepoints != tcs_img.shape[0]:
        #    raise AttributeError('Shape mismatch between list of subjects {} '
        #                         'and timecourses data shape {}.'.format(n_timepoints,
        #                                                                 tcs_img.shape))

    def _load_timecourses(self):
        """ Return the timecourses image file and checks if the shape is correct."""
        # load the timecourses file
        tcsf = fetch_one_file(self.ica_dir, self._tcs_fname)
        tcs = niimg.load_img(tcsf).get_data()
        return tcs

    def _fetch_components_file(self):
        return fetch_one_file(self.ica_dir, self._comps_fname)

    # def _calculate_goodness_of_fit(self):
    #     """ Return the goodness-of-fit values from a GIFT result."""
    #     pass
    #     #TODO
    #
    # def goodness_of_fit_df(self, group_labels_file, subjid_pat=r'(?P<patid>[a-z]{2}_[0-9]{6})'):
    #     """ Return a pandas.DataFrame ready for an excel file with:
    #     - the subject IDs taken from the file paths inside the *Subject.mat file,
    #     - the groups taken from `group_labels_file`, and
    #     - the goodness-of-fit measures for each subject and independent component.
    #
    #     Parameters
    #     ----------
    #     group_labels_file: str
    #         A CSV file with two columns: "subject_id" and "group".
    #         The subject_ids must be in the paths contained in the Subject.mat
    #         file and match the `subjid_pat` argument.
    #
    #     subjid_pat: regext str
    #         A search regex pattern that returns one group element that
    #         contains the subject id.
    #         This will be used to *search* for subject_id in the file paths
    #         contained in the Subjects.mat file.
    #
    #     Returns
    #     -------
    #     gof_df: pandas.DataFrame
    #     """
    #     # make sure file exists
    #     if not op.exists(group_labels_file):
    #         raise FileNotFoundError('The file {} has not been found.'.format(group_labels_file))
    #
    #     # make sure this object has been .fit()
    #     self._update()
    #
    #     # read the groups file
    #     groups = self._parse_groups_file(group_labels_file=group_labels_file)
    #     patids = self._get_subject_ids(subjid_pat=subjid_pat)
    #
    #     # calculate the goodness of fit
    #     gofs = self._calculate_goodness_of_fit(patids)
    #
    #     # build the goodness-of-fit table
    #     df = build_raw_loadings_table(gofs, patids)
    #     df = add_groups_to_loadings_table(df, groups)
    #
    #     return df


class SBMICAResultsPlotter(MIALABICAResultsPlotter):
    """ Use nilearn to plot results from the MIALAB SBM tool, given the ICA result folder path.

    Parameters
    ----------
    ica_result_dir: str
        Path to the ICA output folder or the ICA components volume file.
    """
    _tcs_fname  = '*_loading_coeff_.nii'

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


def ica_loadings_sheet(plotter, labels_file, mask=None, bg_img=None, zscore=2.,
                       subjid_pat=r'(?P<patid>[a-z]{2}_[0-9]{6})'):
    """ Save the Excel loadings files in the `ica_out_dir`.
    One file is `subject_loadings.xls` which has the loadings as is, with the subjects IDs and group.
    The other file is `subject_group_loadings.xls` which has the loading signs changed according to
    the average correlation value of the "main" region of each of the IC spatial maps.

    Parameters
    ----------
    plotter: MIALABICAResultsPlotter, or any derivative
        The GIFT ICA analysis plotter object.
        You must call the plotter fit function before, otherwise this function will call fit
        with `mask`, `zscore` and mode='+-'.

    labels_file: str
        A CSV file with two columns: "subject_id" and "group".
        The subject_ids must be in the paths contained in the Subject.mat
        file and match the `subjid_pat` argument.

    mask: str
        Path to a mask file to select only brain area from the IC spatial map.

    bg_img: str
        A background image for the blob plots check report, to verify that the blobs
        taken into account for the loadings signs are correct.

    zscore: float
        Value to threshold the IC spatial maps to obtain the IC spatial map "main" region.

    subjid_pat: regext str
        A search regex pattern that returns one group element that
        contains the subject id.
        This will be used to search for subject_id in the file paths
        contained in the Subjects.mat file.
    """
    rawloadings_filename   = 'subject_loadings.xls'
    grouploadings_filename = 'subject_weighted_loadings.xls'
    check_blob_plot        = 'check_sign_blobs.png'

    if not plotter.is_fit():
        plotter.fit(mask_file=mask, mode='+-', zscore=zscore)

    output_dir = plotter.ica_dir

    # generate and save the simple loadings sheet
    sdf = plotter.simple_loadings_df(group_labels_file=labels_file, subjid_pat=subjid_pat)
    sdf.to_excel(op.join(output_dir, rawloadings_filename))

    # generate and save the group-processed loadings sheet
    pdf = plotter.weighted_loadings_df(group_labels_file=labels_file, subjid_pat=subjid_pat)
    pdf.to_excel(op.join(output_dir, grouploadings_filename))

    # plot blobs over IC maps for checking
    check_blob_plot = op.join(output_dir, check_blob_plot)
    plotter.plot_icmaps_and_blobs(check_blob_plot, bg_img=bg_img)