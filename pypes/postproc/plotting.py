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


def plot_ica_results(ica_result, application='nilearn', mask_file='', zscore=0, **kwargs):
    """ Use nilearn through pypes to plot results from CanICA and DictLearning, given the ICA result folder path.
    Parameters
    ----------
    ica_result: str
        Path to the ICA output folder or the ICA components volume file.

    application: str
        Choicese: ('nilearn', 'gift')

    mask_file: str
        Path to the brain mask file to be used for thresholding.

    thr: int
        Value of the Z-score thresholding.
    """
    import os.path as op
    from   glob import glob

    import nibabel as nib
    from   boyle.nifti.utils import filter_icc
    from   pypes.nilearn.plot import (plot_all_components,
                                      plot_canica_components)

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
        from nilearn.image import iter_img
        mask = nib.load(mask_file)
        icc_imgs = [filter_icc(icc, zscore, True, mask) for icc in list(iter_img(icc_file))]
    else:
        icc_imgs = icc_file

    # specify the file paths
    all_icc_plot_f = op.join(ica_dir, 'all_components.pdf')
    iccs_plot_f    = op.join(ica_dir, 'canica_components.pdf')

    # make the plots
    fig1 = plot_canica_components(icc_imgs, **kwargs)
    fig1.savefig(iccs_plot_f)

    fig2 = plot_all_components(icc_imgs, **kwargs)
    fig2.savefig(all_icc_plot_f)

    return all_icc_plot_f, iccs_plot_f