# -*- coding: utf-8 -*-
"""
Helper functions to plot results.
"""
import re
import os.path as op
from   glob import glob

import scipy.io           as sio
import pandas             as pd
import nilearn.image      as niimg
from   boyle.nifti.utils  import filter_icc
from   nilearn.image      import iter_img
from   pypes.nilearn.plot import (plot_all_components,
                                  plot_canica_components,
                                  plot_multi_slices,
                                  plot_overlays)

from ..utils.files import fetch_one_file
from .ica_loadings import (filter_ics,
                           get_largest_blobs,
                           build_raw_loadings_table,
                           add_groups_to_loadings_table,
                           )


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


def plot_ica_results(ica_result, application, mask_file='', zscore=2., mode='+', bg_img=''):
    """  Use nilearn to plot results from different ICA analysis tools, given the ICA result folder path."""
    if application == 'gift':
        plotter = GIFTICAResultsPlotter(ica_result)
    elif application == 'canica':
        plotter = CanICAResultsPlotter(ica_result)
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

        icc_files = glob(op.join(self.ica_dir, self._comps_fname))
        if len(icc_files) != 1:
            raise IOError('Expected 1 ICC file, found {}: {}.'.format(len(icc_files), icc_files))
        return icc_files[0]

    def _filter_ic_imgs(self, ic_file):
        # filter the ICC if mask and threshold are set
        if self.mask_file and self.zscore > 0:
            mask = niimg.load_img(self.mask_file)
            icc_imgs = [filter_icc(icc, mask=mask, thr=self.zscore, zscore=True, mode=self.mode)
                        for icc in iter_img(ic_file)]
        else:
            icc_imgs = [niimg.load_img(ic) for ic in ic_file]
        return icc_imgs

    def _update(self, force=False):
        if not hasattr(self, '_icc_imgs'):
            raise RuntimeError('You should call `fit()` before using this function.')

        if self._icc_imgs is not None and not force:
            return

        ic_file = self._fetch_components_file()
        self._icc_imgs = filter_ics(ic_file, mask=self.mask_file, zscore=self.zscore)

    def plot_icmaps(self, outtype='pdf', **kwargs):
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
        icc_multi_slice = op.join(self.ica_dir,  'ic_map_{}_zscore_{}.{}')

        # make the plots
        fig1 = plot_canica_components(self._icc_imgs, **kwargs)
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
                                     title="IC map {} (z-score {})".format(i+1, self.zscore),
                                     title_fontsize=32,
                                     plot_func=None,
                                     **kwargs)
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


class GIFTICAResultsPlotter(ICAResultsPlotter):
    """ Use nilearn to plot results from CanICA and DictLearing tools from nilearn, given the ICA result folder path.

    Parameters
    ----------
    ica_result_dir: str
        Path to the ICA output folder or the ICA components volume file.
    """
    # define the input file patterns
    _comps_fname    = '*_component_ica*.nii'
    _loadings_fname = '*loading_coeff*nii'
    _subjects_fname = '*Subject.mat'

    def simple_loadings_sheet(self, group_labels_file, subjid_pat=r'(?P<patid>[a-z]{2}_[0-9]{6})'):
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
        # # the output file paths
        # rawloadings_filename   = 'subject_loadings.xls'
        # grouploadings_filename = 'subject_group_loadings.xls'
        # check_blob_outdir      = 'check_icmap_blob'

        # make sure this object has been .fit()
        self._update()

        # read the groups file
        groups = None
        if group_labels_file and op.exists(group_labels_file):
            groups = pd.read_csv(group_labels_file)
            if not 'subject_id' in groups.columns or \
               not 'group' in groups.columns:
                raise AttributeError("Please add columns names 'subject_id' and 'group' to the"
                                     "group labels file.")

        # load the .mat file with subjects lists, mainly to get the order in which subjects are
        # introduced in loads
        subjsf   = glob(op.join(self.ica_dir, self._subjects_fname))[0]
        in_files = list(sio.loadmat(subjsf)['files'][0][0][0])

        # process the paths to extract the patient IDs given by `patid_re`
        in_files = [f.split(',')[0] for f in in_files]
        patid_re = re.compile(subjid_pat, re.I)
        patids = [patid_re.search(f).group() for f in in_files]

        # load the loadings file
        loadf = fetch_one_file(self.ica_dir, self._loadings_fname)
        loads = niimg.load_img(loadf).get_data()

        # check shapes
        if len(patids) != loads.shape[0]:
            raise AttributeError('Shape mismatch between list of subjects {} '
                                 'and loadings data shape {}.'.format(len(patids), loads.shape))

        # build the raw loadings table
        df = build_raw_loadings_table(loads, patids)
        df = add_groups_to_loadings_table(df, groups)
        #df.to_excel(op.join(ica_out_dir, rawloadings_filename))

        return df

    def _get_icmaps_blobs(self):
        """ Return the thresholded ic_maps and the blobs in two generators."""
        # make sure this object has been .fit()
        self._update()

        # load the loadings file
        loadf = fetch_one_file(self.ica_dir, self._loadings_fname)
        loads = niimg.load_img(loadf).get_data()

        # step 2: calculate the signs of the components for the groups and save everything in another spreadsheet
        # load the components image file
        compsf = fetch_one_file(self.ica_dir, self._comps_fname)
        comps_img = niimg.load_img(compsf)

        n_ics = loads.shape[1]
        if comps_img.get_data().shape[-1] != n_ics:
            raise AttributeError('Shape mismatch between loadings matrix {} '
                                 'and components data shape {}.'.format(loads.shape, comps_img.get_data().shape))

        return get_largest_blobs(self._icc_imgs)

    def group_loadings_sheet(self, group_labels_file, subjid_pat=r'(?P<patid>[a-z]{2}_[0-9]{6})'):
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
        # make sure this object has been .fit()
        self._update()

        # let's first pick the simple version of the loadings
        df = self.simple_loadings_sheet(group_labels_file, subjid_pat=subjid_pat)
        # group the df by group
        grouped = df.groupby('group')

        # get the values of the largest blobs in the filtered IC maps
        blobs = self._get_icmaps_blobs()

        # calculate the avg per blob
        blob_avgs = [b.get_data().mean() for b in blobs]

        # check if the dimensions match
        group_averages = grouped.mean() #.aggregate(np.mean)
        n_ics = len(blob_avgs)
        if len(group_averages.columns) != n_ics:
            raise AttributeError('The number of IC maps and loadings are not equal, '
                                 'got {} and {}.'.format(n_ics, len(group_averages.columns)))

        # multiply the avg blob value with the group loadings
        sign_df = group_averages * blob_avgs

        # multiply the whole group loadings by the IC blob sign
        for name, group in grouped:
            group[list(range(1, n_ics+1))] = [group[list(range(1, n_ics+1))] * sign_df.ix[name]]

        # put the dataframe back into one whole ungrouped dataframe
        return pd.concat((group for name, group in group.items()))

    def plot_icmaps_and_blobs(self, outdir, outtype='pdf', bg_img=None, **kwargs):
        """ Plot the IC maps with the largest blobs overlaid."""
        # make sure this object has been .fit()
        self._update()
        # TODO
        #fig = plot_overlays(stat_imgs, contour_imgs, bg_img, figsize=(2.5, 3), **kwargs):
        # save fig

        # close fig