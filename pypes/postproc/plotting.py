# -*- coding: utf-8 -*-
"""
Helper functions to plot results.
"""


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


def plot_ica_results(ica_result, application='nilearn', mask_file='', mode='+', zscore=0, **kwargs):
    """ Use nilearn through pypes to plot results from CanICA and DictLearning, given the ICA result folder path.
    Parameters
    ----------
    ica_result: str
        Path to the ICA output folder or the ICA components volume file.

    application: str
        Choicese: ('nilearn', 'gift')

    mask_file: str
        Path to the brain mask file to be used for thresholding.

    mode: str
        Choices: '+' for positive threshold,
                 '+-' for positive and negative threshold and
                 '-' for negative threshold.

    zscore: int or float
        Value of the Z-score thresholding.

    Returns
    -------
    all_icc_plot_f: str

    iccs_plot_f: str

    sliced_ic_plots: list of str
    """
    import os.path as op
    from   glob import glob

    import nibabel as nib
    from   boyle.nifti.utils  import filter_icc
    from   nilearn.image      import iter_img
    from   pypes.nilearn.plot import (plot_all_components,
                                      plot_canica_components,
                                      plot_multi_slices)

    def find_ica_components_file(ica_dir, application='nilearn'):
        base_dir = op.expanduser(ica_dir)

        app_choices = ('nilearn', 'gift')
        if application not in app_choices:
            raise ValueError('Unexpected value for `application` {}, should be any of ({}).'.format(application,
                                                                                                    app_choices))

        if application == 'nilearn':
            fname = 'canica_resting_state.nii.gz'
        elif application == 'gift':
            fname = '*_component_ica*.nii'

        icc_files = glob(op.join(base_dir, fname))

        if len(icc_files) != 1:
            raise IOError('Expected 1 ICC file, found {}: {}.'.format(len(icc_files), icc_files))

        return icc_files[0]

    if not op.exists(ica_result):
        raise IOError('Expected an existing file or folder, but could not find {}.'.format(ica_result))

    if op.isdir(ica_result):
        ica_dir  = ica_result
        icc_file = find_ica_components_file(ica_result, application=application)
    else:
        icc_file = ica_result
        ica_dir  = op.dirname(icc_file)

    # filter the ICC if mask and threshold are set
    if mask_file and zscore > 0:
        mask = nib.load(mask_file)
        icc_imgs = [filter_icc(icc, mask=mask, thr=zscore, zscore=True, mode=mode) for icc in list(iter_img(icc_file))]
    else:
        icc_imgs = icc_file

    # specify the file paths
    all_icc_plot_f  = op.join(ica_dir, 'all_components_zscore_{}.pdf'.format(zscore))
    iccs_plot_f     = op.join(ica_dir, 'ic_components_zscore_{}.pdf'.format(zscore))
    icc_multi_slice = op.join(ica_dir, 'ic_map_{}_zscore_{}.pdf')

    # make the plots
    fig1 = plot_canica_components(icc_imgs, **kwargs)
    fig1.savefig(iccs_plot_f, facecolor=fig1.get_facecolor(), edgecolor='none')

    fig2 = plot_all_components(icc_imgs, **kwargs)
    fig2.savefig(all_icc_plot_f, facecolor=fig2.get_facecolor(), edgecolor='none')

    # make the multi sliced IC plots
    sliced_ic_plots = []
    for i, img in enumerate(iter_img(icc_imgs)):
        fig3 = plot_multi_slices(img,
                                 cut_dir="z",
                                 n_cuts=24,
                                 n_cols=4,
                                 title="IC map {} (z-score {})".format(i+1, zscore),
                                 title_fontsize=32,
                                 plot_func=None,
                                 **kwargs)
        out_f = icc_multi_slice.format(i+1, zscore)
        fig3.savefig(out_f, facecolor=fig3.get_facecolor(), edgecolor='none')
        sliced_ic_plots.append(out_f)

    return all_icc_plot_f, iccs_plot_f, sliced_ic_plots
