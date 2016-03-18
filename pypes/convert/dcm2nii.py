# -*- coding: utf-8 -*-
"""
DICOM to Nifti converter node based on dcm2nii.
"""
import os.path as op

from nipype.interfaces.base    import traits
from nipype.interfaces.dcm2nii import Dcm2nii

from .._utils import format_pair_list
from ..utils  import (spm_tpm_priors_path,
                      extend_trait_list,
                      find_wf_node,
                      get_datasink,
                      get_input_node,
                      remove_ext,
                      get_input_file_name)


def dcm2nii_converter(source_names=traits.Undefined):
    """ Create a dcm2nii workflow.

    It does:
    - dcm2nii

    Nipype Inputs
    -------------
    source_names: traits.File
        path to the DICOM files or series folder.

    Returns
    -------
    wf: nipype Workflow

    Note
    ----
    For more info: http://www.mit.edu/~satra/nipype-nightly/interfaces/generated/nipype.interfaces.dcm2nii.html
    """
    dcm2nii = Dcm2nii()

    dcm2nii.inputs.gzip_output     = True
    dcm2nii.inputs.output_dir      = '.'
    dcm2nii.inputs.terminal_output = 'file'
    dcm2nii.inputs.source_names    = source_names

    return dcm2nii


def attach_dcm2nii(main_wf, wf_name="dcm2nii_convert"):
    """ Attach the SPM12 anatomical MRI pre-processing workflow to the `main_wf`.

    Parameters
    ----------
    main_wf: nipype Workflow

    wf_name: str
        Name of the preprocessing workflow

    Nipype Inputs for `main_wf`
    ---------------------------
    Note: The `main_wf` workflow is expected to have an `input_files` and a `datasink` nodes.

    input_files.select.anat: input node

    datasink: nipype Node

    Returns
    -------
    main_wf: nipype Workflow
    """
    in_files = get_input_node(main_wf)
    datasink = get_datasink  (main_wf)

    # The base name of the 'dicom' file sources for the substitutions
    anat_fbasename = remove_ext(op.basename(get_input_file_name(in_files, 'anat')))

    #
    # # dataSink output substitutions
    # regexp_subst = [
    #                  (r"/{anat}_.*corrected_seg8.mat$", "/{anat}_to_mni_affine.mat"),
    #                  (r"/m{anat}.*_corrected.nii$",     "/{anat}_biascorrected.nii"),
    #                  (r"/wm{anat}.*_corrected.nii$",    "/{anat}_mni.nii"),
    #                  (r"/y_{anat}.*nii$",               "/{anat}_to_mni_field.nii"),
    #                  (r"/iy_{anat}.*nii$",              "/{anat}_to_mni_inv_field.nii"),
    #                  (r"/mwc1{anat}.*nii$",             "/{anat}_gm_mod_mni.nii"),
    #                  (r"/mwc2{anat}.*nii$",             "/{anat}_wm_mod_mni.nii"),
    #                  (r"/mwc3{anat}.*nii$",             "/{anat}_csf_mod_mni.nii"),
    #                  (r"/mwc4{anat}.*nii$",             "/{anat}_nobrain_mod_mni.nii"),
    #                  (r"/c1{anat}.*nii$",               "/{anat}_gm.nii"),
    #                  (r"/c2{anat}.*nii$",               "/{anat}_wm.nii"),
    #                  (r"/c3{anat}.*nii$",               "/{anat}_csf.nii"),
    #                  (r"/c4{anat}.*nii$",               "/{anat}_nobrain.nii"),
    #                  (r"/c5{anat}.*nii$",               "/{anat}_nobrain_mask.nii"),
    #                ]
    # regexp_subst = format_pair_list(regexp_subst, anat=anat_fbasename)
    # datasink.inputs.regexp_substitutions = extend_trait_list(datasink.inputs.regexp_substitutions,
    #                                                          regexp_subst)
    #
    # # input and output anat workflow to main workflow connections
    # main_wf.connect([(in_files, t1_wf,    [("anat",                                  "bias_correction.input_image")]),
    #                  (t1_wf,    datasink, [("warp_anat.normalized_files",            "anat.@mni")],),
    #                  (t1_wf,    datasink, [("new_segment.modulated_class_images",    "anat.tissues.@warped"),
    #                                        ("new_segment.native_class_images",       "anat.tissues.@native"),
    #                                        ("new_segment.transformation_mat",        "anat.transform.@linear"),
    #                                        ("new_segment.forward_deformation_field", "anat.transform.@forward"),
    #                                        ("new_segment.inverse_deformation_field", "anat.transform.@inverse"),
    #                                        ("new_segment.bias_corrected_images",     "anat.@biascor"),
    #                                       ]),
    #                 ])

    return main_wf
