# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:


import nipype.interfaces.spm     as spm
import nipype.pipeline.engine    as pe
from   nipype.algorithms.misc    import Gunzip, TSNR
from   nipype.interfaces.utility import Function, Select, Split, Merge, IdentityInterface
from   nipype.interfaces.nipy.preprocess import Trim
from   nipype.interfaces import fsl

from   ..preproc import (spm_apply_deformations,
                         auto_nipy_slicetime,
                         auto_spm_slicetime,
                         nipy_motion_correction,
                         )

from   .._utils import format_pair_list
from   ..utils import (setup_node,
                       remove_ext,
                       extend_trait_list,
                       get_input_node,
                       get_datasink,
                       get_input_file_name)


def extract_noise_components(realigned_file, noise_mask_file, num_components):
    """Derive components most reflective of physiological noise.
        The noise removal is based on Behzadi et al. (2007)
    """
    import os
    from nibabel import load
    import numpy as np
    import scipy as sp
    imgseries = load(realigned_file)
    components = None
    mask = load(noise_mask_file).get_data()
    voxel_timecourses = imgseries.get_data()[mask > 0]
    voxel_timecourses[np.isnan(np.sum(voxel_timecourses, axis=1)), :] = 0
    # remove mean and normalize by variance
    # voxel_timecourses.shape == [nvoxels, time]
    X = voxel_timecourses.T

    stdX = np.std(X, axis=0)
    stdX[stdX == 0] = 1.
    stdX[np.isnan(stdX)] = 1.
    stdX[np.isinf(stdX)] = 1.
    X = (X - np.mean(X, axis=0)) / stdX
    u, _, _ = sp.linalg.svd(X, full_matrices=False)
    if components is None:
        components = u[:, :num_components]
    else:
        components = np.hstack((components, u[:, :num_components]))
    components_file = os.path.join(os.getcwd(), 'noise_components.txt')
    np.savetxt(components_file, components, fmt="%.10f")
    return components_file


def noise_removal():
    """ Create a noise removal node.

    Nipype Inputs
    -------------


    Nipype Outputs
    --------------


    Returns
    -------
    rmnoise_node: nipype.Node

    """
    rmnoise_input = setup_node(IdentityInterface(
                               fields=["trimmed_rest",],
                               mandatory_inputs=True),
                               name="rmnoise_input")


    tsnr      = setup_node(TSNR(regress_poly=2),
                           name='tsnr')
    getthresh = setup_node(interface=fsl.ImageStats(op_string='-p 98'),
                           name='getthresh')
    threshold_stddev = pe.Node(fsl.Threshold(),
                               name='threshold')
    compcor = pe.Node(Function(input_names=['realigned_file',
                                            'noise_mask_file',
                                            'num_components'],
                                output_names=['noise_components'],
                                function=extract_noise_components),
                      name='compcorr')
    rm_noise = pe.Node(fsl.FilterRegressor(filter_all=True),
                       name='remove_noise')
