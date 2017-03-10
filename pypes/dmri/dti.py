# -*- coding: utf-8 -*-
"""
Diffusion Tensor Image preprocessing and tensor model fitting workflow.
"""
import os.path as op

from .artifacts import attach_dti_artifact_correction
from .coregister import spm_anat_to_diff_coregistration
from .._utils  import format_pair_list
from ..config  import check_atlas_file
from ..utils   import (get_datasink,
                       get_input_node,
                       get_interface_node,
                       remove_ext,
                       extend_trait_list,
                       get_input_file_name,
                       extension_duplicates,
                       )


def attach_spm_fsl_dti_preprocessing(main_wf, wf_name="spm_fsl_dti_preprocessing"):
    """ Attach a set of pipelines to the `main_wf` for Diffusion MR (`diff`) image processing.
    - dti_artifact_correction
    - spm_anat_to_diff_coregistration
    - dti_tensor_fitting

    Parameters
    ----------
    main_wf: nipype Workflow

    wf_name: str
        Name of the preprocessing workflow

    params: dict with parameter values
        atlas_file: str
            Path to the anatomical atlas to be transformed to diffusion MRI space.

    Nipype Inputs for `main_wf`
    ---------------------------
    Note: The `main_wf` workflow is expected to have an `input_files` and a `datasink` nodes.

    input_files.select.diff: input node

    datasink: nipype Node

    Returns
    -------
    main_wf: nipype Workflow
    """
    in_files = get_input_node(main_wf)
    datasink = get_datasink  (main_wf)
    anat_output = get_interface_node(main_wf, 'anat_output')

    # attach the artifact detection and correction pipeline
    main_wf = attach_dti_artifact_correction(main_wf)
    dti_art_output = get_interface_node(main_wf, 'dti_art_output')

    # The workflow boxes
    coreg_dti_wf = spm_anat_to_diff_coregistration(wf_name=wf_name)

    # dataSink output substitutions
    ## The base name of the 'diff' file for the substitutions
    diff_fbasename = remove_ext(op.basename(get_input_file_name(in_files, 'diff')))
    anat_fbasename = remove_ext(op.basename(get_input_file_name(in_files, 'anat')))

    regexp_subst = [
                    (r"/brain_mask_{diff}_space\.nii$", "/brain_mask.nii"),
                    (r"/eddy_corrected\.nii$",          "/{diff}_eddycor.nii"),
                    (r"/rc1anat_hc_corrected\.nii$",    "/gm_diff.nii"),
                    (r"/rc2anat_hc_corrected\.nii$",    "/wm_diff.nii"),
                    (r"/rc3anat_hc_corrected\.nii$",    "/csf_diff.nii"),
                    (r"/rmanat_hc_corrected\.nii$",     "/{anat}_diff.nii"),
                   ]
    regexp_subst = format_pair_list(regexp_subst, diff=diff_fbasename,
                                                  anat=anat_fbasename)

    # prepare substitution for atlas_file, if any
    do_atlas, atlas_file = check_atlas_file()
    if do_atlas:
        atlas_basename = remove_ext(op.basename(atlas_file))
        regexp_subst.extend([
                             (r"/[\w]*{atlas}.*\.nii$", "/{atlas}_{diff}_space.nii"),
                            ])
        regexp_subst = format_pair_list(regexp_subst, atlas=atlas_basename,
                                                      diff=diff_fbasename)


    regexp_subst += extension_duplicates(regexp_subst)
    datasink.inputs.regexp_substitutions = extend_trait_list(datasink.inputs.regexp_substitutions,
                                                             regexp_subst)

    # input and output diffusion MRI workflow to main workflow connections
    main_wf.connect([
                     (dti_art_output, coreg_dti_wf, [("avg_b0",         "dti_co_input.avg_b0"),]),
                     (anat_output,    coreg_dti_wf, [("tissues_native", "dti_co_input.tissues"),
                                                     ("anat_biascorr",  "dti_co_input.anat")
                                                    ]),
                     (coreg_dti_wf, datasink, [("dti_co_output.anat_diff",       "diff.@anat_diff"),
                                               ("dti_co_output.tissues_diff",    "diff.tissues.@tissues_diff"),
                                               ("dti_co_output.brain_mask_diff", "diff.@brain_mask"),
                                              ]),
                    ])

    if do_atlas:
            main_wf.connect([(anat_output,  coreg_dti_wf, [("atlas_anat", "dti_co_input.atlas_anat")]),
                             (coreg_dti_wf, datasink,     [("dti_co_output.atlas_diff", "diff.@atlas")]),
                            ])

    return main_wf


