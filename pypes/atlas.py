# -*- coding: utf-8 -*-
"""
Nipype workflows to process brain atlases.
"""
import os.path as op

import nipype.interfaces.spm     as spm
import nipype.pipeline.engine    as pe
from   nipype.algorithms.misc    import Gunzip
from   nipype.interfaces.ants    import N4BiasFieldCorrection
from   nipype.interfaces.base    import traits
from   nipype.interfaces.io      import DataSink, SelectFiles


def spm_apply_inverse_transformation(in_imgs=traits.Undefined, trans_field=traits.Undefined):


def spm_apply_deformations(in_imgs=traits.Undefined, trans_field=traits.Undefined):
    """Return a Normalize12 interface object.
    SPM12's new Normalise routine for warping an image to a template.
    For more info:
    http://www.mit.edu/~satra/nipype-nightly/interfaces/generated/nipype.interfaces.spm.preprocess.html#normalize12

    Parameters
    ----------
    in_imgs: iterable of str

    trans_field:
        file y_*.nii containing 3 deformation fields for the deformation in
        x, y and z dimension

    Ouputs
    ------
    deformation_field: (a list of items which are an existing file name)
        NIfTI file containing 3 deformation fields for the deformation in x,
        y and z dimension

    normalized_files: (a list of items which are an existing file name)
        Normalized other files

    normalized_image: (a list of items which are an existing file name)
        Normalized file that needed to be aligned

    """
    norm12 = spm.Normalize12(jobtype='write')
    norm12.inputs.deformation_file  = trans_field
    norm12.inputs.image_to_align    = in_imgs
    norm12.inputs.write_voxel_sizes = [1, 1, 1]

    #norm12.run()
    return norm12
