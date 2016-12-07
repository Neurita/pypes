# -*- coding: utf-8 -*-
"""
Nipype utilities for processing anatomical MRI.
"""
import nipype.interfaces.spm as spm
from   nipype.interfaces.ants import N4BiasFieldCorrection
from   nipype.interfaces.base import traits

from   ..utils import spm_tpm_priors_path


def biasfield_correct(anat_filepath=traits.Undefined):
    """ Inhomogeneity correction.
    ANTS N4BiasFieldCorrection interface.

    Parameters
    ----------
    anat_filepath: str
        Path to the anatomical file path

    Returns
    -------
    seg: N4BiasFieldCorrection interface
    """
    n4 = N4BiasFieldCorrection()
    n4.inputs.dimension = 3
    n4.inputs.bspline_fitting_distance = 300
    n4.inputs.shrink_factor = 3
    n4.inputs.n_iterations = [50, 50, 30, 20]
    n4.inputs.convergence_threshold = 1e-6
    #n4.inputs.bspline_order = 5
    n4.inputs.save_bias = True

    n4.inputs.input_image = anat_filepath

    return n4


def spm_segment(anat_filepath=traits.Undefined, priors_path=None):
    """ SPM12 New Segment interface.

    Parameters
    ----------
    anat_filepath: str
        Path to the anatomical file path

    priors_path: str
        Path to the tissue probability maps file

    Returns
    -------
    seg: NewSegment interface
    """
    if priors_path is None:
        priors_path = spm_tpm_priors_path()

    seg = spm.NewSegment()

    tissue1 = ((priors_path, 1), 1, (True,  True),   (True,  True))
    tissue2 = ((priors_path, 2), 1, (True,  True),   (True,  True))
    tissue3 = ((priors_path, 3), 2, (True,  True),   (True,  True))
    tissue4 = ((priors_path, 4), 3, (True,  True),   (True,  True))
    tissue5 = ((priors_path, 5), 4, (True,  False),  (False, False))
    tissue6 = ((priors_path, 6), 2, (False, False),  (False, False))
    seg.inputs.tissues = [tissue1, tissue2, tissue3, tissue4, tissue5, tissue6]
    seg.inputs.channel_info = (0.0001, 60, (True, True))
    #seg.inputs.warping_regularization = [0, 0.001, 0.5, 0.05, 0.2]
    seg.inputs.write_deformation_fields = [True, True]

    seg.inputs.channel_files = anat_filepath

    #seg.run()
    return seg
