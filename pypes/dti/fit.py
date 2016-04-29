# -*- coding: utf-8 -*-
"""
Nipype workflows to preprocess diffusion MRI.
"""
import os.path as op

import nipype.pipeline.engine    as pe
from   nipype.interfaces.fsl     import ExtractROI, Eddy, MultiImageMaths
from   nipype.interfaces.utility import Function, Select, Split, IdentityInterface
from   nipype.algorithms.misc    import Gunzip
from   nipype.workflows.dmri.fsl.utils import eddy_rotate_bvecs

from   .utils import dti_acquisition_parameters

from   .._utils  import flatten_list, format_pair_list
from   ..preproc import spm_coregister
from   ..config  import setup_node, check_atlas_file
from   ..utils   import (get_datasink,
                         get_input_node,
                         remove_ext,
                         extend_trait_list,
                         get_input_file_name,
                         extension_duplicates,
                         )


def fsl_dti_preprocessing(wf_name="fsl_dti_preproc"):
    """ Run the diffusion MRI pre-processing workflow against the diff files in `data_dir`.

    This estimates an affine transform from anat to diff space, applies it to
    the brain mask and an atlas, and performs eddy-current correction.
    Nipype Inputs
    -------------
    dti_input.diff: traits.File
        path to the diffusion MRI image

    dti_input.bval: traits.File
        path to the bvals file

    dti_input.bvec: traits.File
        path to the bvecs file

    dti_input.tissues: traits.File
        paths to the NewSegment c*.nii output files

    dti_input.anat: traits.File
        path to the high-contrast anatomical image

    dti_input.atlas_anat: traits.File
        Atlas in subject anatomical space.

    Nipype Outputs
    --------------
    dti_output.diff_corrected: traits.File
        Eddy currents corrected DTI image.

    dti_output.bvec_rotated: traits.File
        Rotated bvecs file

    dti_output.brain_mask_diff: traits.File
        Brain mask for diffusion image.

    dti_output.atlas_diff: traits.File
        Atlas image warped to diffusion space.
        If the `atlas_file` option is an existing file and `normalize_atlas` is True.

    Nipype Workflow Dependencies
    ----------------------------
    This workflow depends on:
    - spm_anat_preproc

    Returns
    -------
    wf: nipype Workflow
    """

    dti_input    = setup_node(IdentityInterface(
        fields=["diff", "bval", "bvec", "tissues", "anat", "atlas_anat"],
        mandatory_inputs=True),
        name="dti_input")

    write_acqp   = setup_node(Function(function=dti_acquisition_parameters,
                                       input_names=["in_file"],
                                       output_names=["out_acqp", "out_index"],
        ),                    name="write_acqp")
    extract_b0   = setup_node(ExtractROI(t_min=0, t_size=1),         name="extract_b0")
    gunzip_b0    = setup_node(Gunzip(),                              name="gunzip_b0")
    coreg_b0     = setup_node(spm_coregister(cost_function="mi"),    name="coreg_b0")
    brain_sel    = setup_node(Select(index=[0, 1, 2]),               name="brain_sel")
    coreg_split  = setup_node(Split(splits=[1, 2], squeeze=True),    name="coreg_split")
    brain_merge  = setup_node(MultiImageMaths(),                     name="brain_merge")
    eddy         = setup_node(Eddy(),                                name="eddy")
    rot_bvec     = setup_node(Function(function=eddy_rotate_bvecs,
                                       input_names=["in_bvec", "eddy_params"],
                                       output_names=["out_file"],),
                              name="rot_bvec")
    dti_output   = setup_node(IdentityInterface(fields=["diff_corrected",
                                                        "atlas_diff",
                                                        "bvec_rotated",
                                                        "bval",
                                                        "brain_mask_diff"]),
                              name="dti_output")

    brain_merge.inputs.op_string = "-add '%s' -add '%s' -abs -kernel gauss 2 -dilM -ero -bin"
    brain_merge.inputs.out_file = "brain_mask_diff_space.nii.gz"

    # Create the workflow object
    wf = pe.Workflow(name=wf_name)

    # Connect the nodes
    wf.connect([
                (dti_input,     extract_b0,     [("diff",                "in_file")]),
                (dti_input,     write_acqp,     [("diff",                "in_file")]),
                (dti_input,     eddy,           [("diff",                "in_file")]),
                (dti_input,     eddy,           [("bval",                "in_bval")]),
                (dti_input,     eddy,           [("bvec",                "in_bvec")]),
                (dti_input,     rot_bvec,       [("bvec",                "in_bvec")]),
                (dti_input,     brain_sel,      [("tissues",             "inlist")]),
                (dti_input,     coreg_b0,       [("anat",                "source")]),
                (write_acqp,    eddy,           [("out_acqp",            "in_acqp"),
                                                 ("out_index",           "in_index")]),
                (extract_b0,    gunzip_b0,      [("roi_file",            "in_file")]),
                (gunzip_b0,     coreg_b0,       [("out_file",            "target")]),
                (coreg_b0,      coreg_split,    [("coregistered_files",  "inlist")]),
                (brain_sel,     coreg_b0,       [(("out", flatten_list), "apply_to_files")]),
                (coreg_split,   brain_merge,    [("out1",                "in_file")]),
                (coreg_split,   brain_merge,    [("out2",                "operand_files")]),
                (brain_merge,   eddy,           [("out_file",            "in_mask")]),
                (brain_merge,   dti_output,     [("out_file",            "brain_mask_diff")]),
                (eddy,          dti_output,     [("out_corrected",       "diff_corrected")]),
                (eddy,          rot_bvec,       [("out_parameter",       "eddy_params")]),
                (rot_bvec,      dti_output,     [("out_file",            "bvec_rotated")]),
                (dti_input,     dti_output,     [("bval",                "bval")]),
              ])

    # add more nodes if to perform atlas registration
    do_atlas, _ = check_atlas_file()
    if do_atlas:
        coreg_atlas = setup_node(spm_coregister(cost_function="mi"), name="coreg_atlas")

        # set the registration interpolation to nearest neighbour.
        coreg_atlas.inputs.write_interp = 0
        wf.connect([
                    (dti_input,   coreg_atlas, [("anat",               "source")]),
                    (gunzip_b0,   coreg_atlas, [("out_file",           "target")]),
                    (dti_input,   coreg_atlas, [("atlas_anat",         "apply_to_files")]),
                    (coreg_atlas, dti_output,  [("coregistered_files", "atlas_diff")]),
                  ])

    return wf


def attach_fsl_dti_preprocessing(main_wf, wf_name="fsl_dti_preproc"):
    """ Attach the FSL-based diffusion MRI pre-processing workflow to the `main_wf`.

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
    anat_wf  = main_wf.get_node("spm_anat_preproc")

    # The workflow box
    dti_wf = fsl_dti_preprocessing(wf_name=wf_name)

    # dataSink output substitutions
    ## The base name of the 'diff' file for the substitutions
    diff_fbasename = remove_ext(op.basename(get_input_file_name(in_files, 'diff')))

    regexp_subst = [
                    (r"/brain_mask_{diff}_space\.nii$", "/brain_mask.nii"),
                    (r"/eddy_corrected\.nii$",          "/{diff}_eddycor.nii"),
                   ]
    regexp_subst = format_pair_list(regexp_subst, diff=diff_fbasename)

    # prepare substitution for atlas_file, if any
    do_atlas, atlas_file = check_atlas_file()
    if do_atlas:
        atlas_basename = remove_ext(op.basename(atlas_file))
        regexp_subst.extend([
                             (r"/[\w]*{atlas}\.nii$", "/{atlas}_{diff}_space.nii"),
                            ])
        regexp_subst = format_pair_list(regexp_subst, atlas=atlas_basename,
                                                      diff=diff_fbasename)


    regexp_subst += extension_duplicates(regexp_subst)
    datasink.inputs.regexp_substitutions = extend_trait_list(datasink.inputs.regexp_substitutions,
                                                             regexp_subst)

    # input and output diffusion MRI workflow to main workflow connections
    main_wf.connect([(in_files, dti_wf,   [("diff",                                  "dti_input.diff"),
                                           ("bval",                                  "dti_input.bval"),
                                           ("bvec",                                  "dti_input.bvec")
                                          ]),
                     (anat_wf,  dti_wf,   [("new_segment.native_class_images",       "dti_input.tissues"),
                                           ("new_segment.bias_corrected_images",     "dti_input.anat")
                                          ]),
                     (dti_wf,   datasink, [("dti_output.diff_corrected",             "diff.@eddy_corrected"),
                                           ("dti_output.brain_mask_diff",            "diff.@mask"),
                                           ("dti_output.bvec_rotated",               "diff.@bvec_rotated"),
                                           ("dti_output.bval",                       "diff.@bval")
                                          ]),
                    ])

    if do_atlas:
            main_wf.connect([(anat_wf, dti_wf,   [("anat_output.atlas_warped",  "dti_input.atlas_anat")]),
                             (dti_wf,  datasink, [("dti_output.atlas_diff",     "diff.@atlas")]),
                            ])

    return main_wf
