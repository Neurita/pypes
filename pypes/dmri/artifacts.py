# -*- coding: utf-8 -*-
"""
Nipype workflows to preprocess diffusion MRI.
"""
import os.path as op

import nipype.pipeline.engine    as pe
from   nipype.interfaces.fsl     import Eddy, BET, ExtractROI
import nipype.interfaces.dipy as dipy
from   nipype.interfaces.utility import Function, Select, Split, IdentityInterface
from   nipype.algorithms.misc    import Gunzip
from   nipype.workflows.dmri.fsl.utils import eddy_rotate_bvecs, b0_average, b0_indices
from   nipype.workflows.dmri.fsl import hmc_pipeline

from   .utils import dti_acquisition_parameters, rapidart_dti_artifact_detection

from   .._utils  import format_pair_list
from   ..config  import setup_node, check_atlas_file, get_config_setting
from   ..utils   import (get_datasink,
                         get_input_node,
                         remove_ext,
                         extend_trait_list,
                         get_input_file_name,
                         extension_duplicates,
                         )


def dti_artifact_correction(wf_name="dti_artifact_correction"):
    """ Run the diffusion MRI pre-processing workflow against the diff files in `data_dir`.

    It will resample/regrid the diffusion image to have isometric voxels.
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
    dti_output.diff_corrected: traits.File
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
    out_fields = ["diff_corrected",
                  "bvec_rotated",
                  "brain_mask_diff",
                  "acqp",
                  "index"]

    # input interface
    dti_input = setup_node(IdentityInterface(fields=in_fields, mandatory_inputs=True),
                           name="dti_input")

    # resample
    resample = setup_node(dipy.Resample(), name='dti_resample')

    # head motion correction
    list_b0 = pe.Node(Function(function=b0_indices,
                               input_names=['in_bval'],
                               output_names=['out_idx'],),
                               name='b0_indices')

    ## extract b0s for brain mask
    extract_b0 = pe.Node(ExtractROI(t_size=1),
                         name="extract_b0")

    bet_dwi0 = pe.Node(BET(frac=0.3, mask=True, robust=True),
                       name='bet_dwi_pre')

    pick_first = lambda lst: lst[0]

    hmc = hmc_pipeline()

    # motion artifacts detection
    do_rapidart = get_config_setting("dmri.artifact_detect", True)
    if do_rapidart:
        art = setup_node(rapidart_dti_artifact_detection(), name="detect_artifacts")

    # Eddy
    eddy = setup_node(Eddy(method='jac'), name="eddy")

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
                # resample to iso-voxel
                (dti_input, resample, [("diff", "in_file"),]),

                # read from input file the acquisition parameters
                (dti_input, write_acqp, [("diff", "in_file")]),

                # reference mask for hmc
                (list_b0,  extract_b0, [(("out_idx", pick_first), "t_min"),]),
                (resample, extract_b0, [("out_file",              "in_file")]),

                (extract_b0, bet_dwi0, [("out_file", "in_file")]),

                # head motion correction (hmc)
                (dti_input, list_b0, [("bval",      "in_bval"),]),
                (dti_input, hmc,     [("bval",      "inputnode.in_bval"),
                                      ("bvec",      "inputnode.in_bvec"),
                                     ]),
                (resample,  hmc,     [("out_file",  "inputnode.in_file")]),
                (bet_dwi0,  hmc,     [("mask_file", "inputnode.in_mask")]),
                (list_b0,   hmc,     [(("out_idx", pick_first), "inputnode.ref_num"),]),

                # Eddy
                #(brain_merge,   eddy,           [("out_file",            "in_mask")]),
                (dti_input,     eddy,           [("diff",                "in_file")]),
                (dti_input,     eddy,           [("bval",                "in_bval")]),
                (dti_input,     eddy,           [("bvec",                "in_bvec")]),
                (write_acqp,    eddy,           [("out_acqp",            "in_acqp"),
                                                 ("out_index",           "in_index")]),

                # rotate bvecs
                (dti_input,     rot_bvec,       [("bvec",                "in_bvec")]),
                (eddy,          rot_bvec,       [("out_parameter",       "eddy_params")]),

                # output
                (write_acqp,    dti_output,     [("out_acqp",            "acqp"),
                                                 ("out_index",           "index")]),
                (eddy,          dti_output,     [("out_corrected",       "eddy_corrected")]),
                (rot_bvec,      dti_output,     [("out_file",            "bvec_rotated")]),
                (dti_input,     dti_output,     [("bval",                "bval")]),
              ])

    if do_rapidart:
        wf.connect([
                    # artifact detection
                    (dti_input, art, [("in_file",        "realigned_files"),
                                      ("motion_params",  "realignment_parameters"),
                                      ("brain_mask",     "mask_file"),
                                     ]),

                    # output
                    (art, dti_output, [("displacement_files",   "art_displacement_files"),
                                       ("intensity_files",      "art_intensity_files"),
                                       ("norm_files",           "art_norm_files"),
                                       ("outlier_files",        "art_outlier_files"),
                                       ("plot_files",           "art_plot_files"),
                                       ("statistic_files",      "art_statistic_files"),
                                      ]),
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
