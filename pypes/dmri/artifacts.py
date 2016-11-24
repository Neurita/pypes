# -*- coding: utf-8 -*-
"""
Nipype workflows to preprocess diffusion MRI.
"""
import os.path as op

import nipype.pipeline.engine    as pe
from   nipype.interfaces.fsl     import Eddy, MultiImageMaths
from   nipype.interfaces.utility import Function, Select, Split, IdentityInterface
from   nipype.algorithms.misc    import Gunzip
from   nipype.workflows.dmri.fsl.utils import eddy_rotate_bvecs, b0_average
from   nipype.workflows.dmri.fsl import hmc_pipeline

from   .utils import dti_acquisition_parameters, nlmeans_denoise, rapidart_dti_artifact_detection

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


def dti_artifact_correction(wf_name="fsl_dti_preproc"):
    """ Run the diffusion MRI pre-processing workflow against the diff files in `data_dir`.

    Corrects for head motion correction and Eddy currents.
    Estimates motion outliers and exports motion reports using nipype.algorithms.RapidArt.

    Nipype Inputs
    -------------
    dti_input.diff: traits.File
        path to the diffusion MRI image

    dti_input.bval: traits.File
        path to the bvals file

    dti_input.bvec: traits.File
        path to the bvecs file


    Nipype Outputs
    --------------
    dti_output.corrected: traits.File
        Eddy currents corrected DTI image.

    dti_output.bvec_rotated: traits.File
        Rotated bvecs file

    dti_output.brain_mask_diff: traits.File
        Brain mask for diffusion image.

    dti_output.acpq: traits.File
        Text file with acquisition parameters calculated for Eddy.

    dti_output.index: traits.File
        Text file with acquisition indices calculated for Eddy.

    Returns
    -------
    wf: nipype Workflow
    """
    # specify input and output fields
    in_fields  = ["diff", "bval", "bvec"]
    out_fields = ["corrected",
                  "bvec_rotated",
                  "brain_mask_diff",
                  "acqp",
                  "index"]

    # input interface
    dti_input = setup_node(IdentityInterface(fields=in_fields, mandatory_inputs=True),
                           name="dti_input")

    # head motion correction
    hm = hmc_pipeline()

    # motion artifacts detection
    art = setup_node(rapidart_dti_artifact_detection(), name="detect_artifacts")

    # Eddy
    eddy = setup_node(Eddy(method='jac'), name="eddy")

    ## extract b0s for brain mask
    extract_b0 = setup_node(Function(function=b0_average,
                                     input_names=['in_dwi', 'in_bval'],
                                     output_names=['out_file']),
                            name='extract_b0')

    #extract_b0   = pe.Node   (ExtractROI(t_min=0, t_size=1),      name="extract_b0")
    gunzip_b0    = pe.Node   (Gunzip(),                           name="gunzip_b0")
    coreg_b0     = setup_node(spm_coregister(cost_function="mi"), name="coreg_b0")

    ## acquisition parameters for Eddy
    write_acqp = setup_node(Function(function=dti_acquisition_parameters,
                                     input_names=["in_file"],
                                     output_names=["out_acqp", "out_index"],),
                            name="write_acqp")

    ## rotate b-vecs
    rot_bvec = setup_node(Function(function=eddy_rotate_bvecs,
                                   input_names=["in_bvec", "eddy_params"],
                                   output_names=["out_file"],),
                          name="rot_bvec")

    # output interface
    dti_output = setup_node(IdentityInterface(fields=out_fields),
                            name="dti_output")

    # Create the workflow object
    wf = pe.Workflow(name=wf_name)

    # Connect the nodes
    wf.connect([
                #

                # head motion correction

                # average b0 images
                (dti_input,     extract_b0,     [("diff",                "in_dwi"),
                                                 ("bval",                "in_bval"),
                                                ]),

                (dti_input,     write_acqp,     [("diff",                "in_file")]),

                # artifact detection
                (dti_input,     art,            [("in_file",        "realigned_files"),
                                                 ("motion_params",  "realignment_parameters"),
                                                 ("brain_mask",     "mask_file"),
                                                ]),

                # Co-registration
                (dti_input,     coreg_b0,       [("anat",                "source")]),

                (extract_b0,    gunzip_b0,      [("out_file",            "in_file")]),
                (coreg_b0,      coreg_split,    [("coregistered_files",  "inlist")]),
                (coreg_split,   brain_merge,    [("out1",                "in_file")]),
                (coreg_split,   brain_merge,    [("out2",                "operand_files")]),

                # Eddy
                (brain_merge,   eddy,           [("out_file",            "in_mask")]),
                (dti_input,     eddy,           [("diff",                "in_file")]),
                (dti_input,     eddy,           [("bval",                "in_bval")]),
                (dti_input,     eddy,           [("bvec",                "in_bvec")]),
                (write_acqp,    eddy,           [("out_acqp",            "in_acqp"),
                                                 ("out_index",           "in_index")]),

                # rotate bvecs
                (dti_input,     rot_bvec,       [("bvec",                "in_bvec")]),
                (eddy,          rot_bvec,       [("out_parameter",       "eddy_params")]),

                # nlmeans denoise
                (eddy,          denoise,        [("out_corrected",       "in_file")]),
                (brain_merge,   denoise,        [("out_file",            "mask_file")]),

                # output
                (write_acqp,    dti_output,     [("out_acqp",            "acqp"),
                                                 ("out_index",           "index")]),
                (brain_merge,   dti_output,     [("out_file",            "brain_mask_diff")]),
                (denoise,       dti_output,     [("out_file",            "denoised")]),
                (eddy,          dti_output,     [("out_corrected",       "eddy_corrected")]),
                (denoise,       dti_output,     [("out_file",            "corrected")]),
                (rot_bvec,      dti_output,     [("out_file",            "bvec_rotated")]),
                (dti_input,     dti_output,     [("bval",                "bval")]),
              ])

    # add more nodes if to perform atlas registration
    if do_atlas:
        coreg_atlas = setup_node(spm_coregister(cost_function="mi"), name="coreg_atlas")

        # set the registration interpolation to nearest neighbour.
        coreg_atlas.inputs.write_interp = 0
        wf.connect([
                    (dti_input,   coreg_atlas, [("anat",               "source"),
                                                ("atlas_anat",         "apply_to_files"),
                                               ]),
                    (gunzip_b0,   coreg_atlas, [("out_file",           "target")]),
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
                             (r"/[\w]*{atlas}.*\.nii$", "/{atlas}_{diff}_space.nii"),
                            ])
        regexp_subst = format_pair_list(regexp_subst, atlas=atlas_basename,
                                                      diff=diff_fbasename)


    regexp_subst += extension_duplicates(regexp_subst)
    datasink.inputs.regexp_substitutions = extend_trait_list(datasink.inputs.regexp_substitutions,
                                                             regexp_subst)

    # input and output diffusion MRI workflow to main workflow connections
    main_wf.connect([(in_files, dti_wf,   [("diff",                              "dti_input.diff"),
                                           ("bval",                              "dti_input.bval"),
                                           ("bvec",                              "dti_input.bvec")
                                          ]),
                     (anat_wf,  dti_wf,   [("new_segment.native_class_images",   "dti_input.tissues"),
                                           ("new_segment.bias_corrected_images", "dti_input.anat")
                                          ]),
                     (dti_wf,   datasink, [("dti_output.eddy_corrected",         "diff.@eddy_corrected"),
                                           ("dti_output.denoised",               "diff.@denoised"),
                                           ("dti_output.corrected",              "diff.@corrected"),
                                           ("dti_output.brain_mask_diff",        "diff.@mask"),
                                           ("dti_output.bvec_rotated",           "diff.@bvec_rotated"),
                                           ("dti_output.bval",                   "diff.@bval"),
                                           ("dti_output.acqp",                   "diff.@acqp"),
                                           ("dti_output.index",                  "diff.@index")
                                          ]),
                    ])

    if do_atlas:
            main_wf.connect([(anat_wf, dti_wf,   [("anat_output.atlas_anat",   "dti_input.atlas_anat")]),
                             (dti_wf,  datasink, [("dti_output.atlas_diff",    "diff.@atlas")]),
                            ])

    return main_wf
