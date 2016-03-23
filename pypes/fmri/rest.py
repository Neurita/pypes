# -*- coding: utf-8 -*-
"""
Nipype workflows to process anatomical MRI.
"""
import os.path as op

import nipype.interfaces.spm     as spm
import nipype.pipeline.engine    as pe
from   nipype.algorithms.misc    import Gunzip, TSNR
from   nipype.interfaces.utility import Function, Select, Split, Merge, IdentityInterface
from   nipype.interfaces         import fsl
from   nipype.interfaces.io      import add_traits
from   nipype.interfaces.nipy.preprocess import Trim, ComputeMask

from   ..preproc import (spm_apply_deformations,
                         auto_nipy_slicetime,
                         auto_spm_slicetime,
                         nipy_motion_correction,
                         extract_noise_components,
                         spm_coregister,
                         rest_noise_filter_wf,
                         bandpass_filter,
                         )


from   .._utils import format_pair_list, flatten_list
from   ..utils import (setup_node,
                       remove_ext,
                       check_atlas_file,
                       extend_trait_list,
                       get_input_node,
                       get_datasink,
                       get_input_file_name,
                       extension_duplicates)


def rest_preprocessing_wf(wf_name="rest_preproc"):
    """ Run the resting-state fMRI pre-processing workflow against the rest files in `data_dir`.

    It does:
    - Trim first 6 volumes of the rs-fMRI file.
    - Slice Timing Correction
    - Atlas in anatomical space coregistration to
    fMRI space (if atlas_file is set in config)

    Nipype Inputs
    -------------
    rest_input.in_files: traits.File
        path to the resting-state image

    rest_input.anat: traits.File
        path to the high-contrast anatomical image

    rest_input.coreg_target: traits.File
        path to the bias corrected anatomical image for coregistration.

    rest_input.tissues: traits.File
        path to the tissue segmentations in anatomical space.

    rest_input.highpass_sigma:traits.Float
        Band pass timeseries filter higher bound in Hz.

    rest_input.lowpass_sigma: traits.Float
        Band pass timeseries filter lower bound in Hz.

    Nipype Outputs
    --------------
    rest_output.motion_corrected: traits.File
        The motion corrected file.

    rest_output.motion_params: traits.File
        The affine transformation file.

    rest_output.time_filtered: traits.File
        The bandpass time filtered fMRI file.

    rest_output.smooth: traits.File
        The isotropic smoothed time filtered image.

    rest_output.epi_brain_mask: traits.File
        An estimated brain mask from mean EPI volume.

    rest_output.tissues_brain_mask: traits.File
        A brain mask calculated from the addition of coregistered
        GM, WM and CSF segmentation volumes from the anatomical
        segmentation.

    rest_output.tissues: traits.File
        The tissues segmentation volume in fMRI space.

    rest_output.anat: traits.File
        The T1w image in fMRI space.

    rest_output.atlas: traits.File
        Atlas image warped to fMRI space.
        If the `atlas_file` option is an existing file and `normalize_atlas` is True.

    Returns
    -------
    wf: nipype Workflow
    """
    # Create the workflow object
    wf = pe.Workflow(name=wf_name)

    # input identities
    rest_input = setup_node(IdentityInterface(fields=["in_files",
                                                      "anat",
                                                      "coreg_target",
                                                      "tissues",
                                                      "atlas_anat",
                                                      "lowpass_freq",
                                                      "highpass_freq",
                                                      ], mandatory_inputs=True),
                            name="rest_input")

    # rs-fMRI preprocessing nodes
    trim    = setup_node(Trim(),
                         name="trim")
    stc_wf  = auto_spm_slicetime()
    realign = setup_node(nipy_motion_correction(), name='realign')

    # Take the mean over time to get a target volume
    meanvol     = setup_node(fsl.MeanImage(), name="meanvol")
    mean_gunzip = setup_node(Gunzip(),        name="mean_gunzip")

    # co-registration nodes
    coreg     = setup_node(spm_coregister(cost_function="mi"), name="coreg")
    brain_sel = setup_node(Select(index=[0, 1, 2]),            name="brain_sel")

    # brain masks
    epi_mask     = setup_node(ComputeMask(), name='epi_mask')
    tissue_mask  = setup_node(fsl.MultiImageMaths(), name='tissue_mask')
    tissue_mask.inputs.op_string = "-add %s -add %s -kernel gauss 2 -dilM -ero -bin"
    gm_select    = setup_node(Select(index=[0]),     name="gm_sel")
    wmcsf_select = setup_node(Select(index=[1, 2]),  name="wmcsf_sel")

    # noise filter
    noise_wf = rest_noise_filter_wf()

    # bandpass filtering
    bandpass = setup_node(Function(input_names=['files', 'lowpass_freq',
                                                'highpass_freq', 'fs'],
                                   output_names=['out_files'],
                                   function=bandpass_filter),
                    name='bandpass_filter')

    # smoothing
    smooth = setup_node(interface=fsl.IsotropicSmooth(fwhm=8), name="fmri_smooth")

    # output identities
    rest_output = setup_node(IdentityInterface(fields=[
                                                       "motion_corrected",
                                                       "motion_params",
                                                       "tissues",
                                                       "anat",
                                                       "time_filtered",
                                                       "smooth",
                                                       "epi_brain_mask",
                                                       "tissues_brain_mask",
                                                       ],
                                               mandatory_inputs=True),
                             name="rest_output")

    # atlas registration?
    do_atlas, atlas_file = check_atlas_file()
    if do_atlas:
        add_traits(rest_input.inputs,  "atlas_anat")
        add_traits(rest_output.inputs, "atlas_rest")
        rest_input.inputs.set(atlas_file=atlas_file)

    # Connect the nodes
    wf.connect([
                # trim
                (rest_input,   trim,         [("in_files", "in_file")]),

                #slice time correction
                (trim,        stc_wf,        [("out_file", "stc_input.in_file")]),

                # motion correction
                (stc_wf,      realign,       [("stc_output.timecorrected_files", "in_file")]),

                # coregistration target
                (realign,     meanvol,       [("out_file", "in_file")]),
                (meanvol,     mean_gunzip,   [("out_file", "in_file")]),
                (mean_gunzip, coreg,         [("out_file", "target")]),

                # epi brain mask
                (mean_gunzip, epi_mask,      [("out_file", "mean_volume")]),

                # coregistration
                (rest_input,  coreg,         [("anat",      "source")]),
                (rest_input,  brain_sel,     [("tissues",   "inlist")]),

                # tissue brain mask
                (coreg,        gm_select,    [("coregistered_files",  "inlist")]),
                (coreg,        wmcsf_select, [("coregistered_files",  "inlist")]),
                (gm_select,    tissue_mask,  [(("out", flatten_list), "in_file")]),
                (wmcsf_select, tissue_mask,  [(("out", flatten_list), "operand_files")]),

                # motion statistics


                # nuisance correction


                # median angle?

                # temporal filtering
                (stc_wf,     bandpass,     [("stc_output.time_repetition", "tr")]),
                (rest_input, bandpass,     [("lowpass_freq",               "lowpass_freq"),
                                            ("highpass_freq",              "highpass_freq"),
                                           ]),
                (realign,    bandpass,     [("out_file",                   "files")]),

                # smoothing
                (bandpass,    smooth,      [("out_files",                   "in_file")]),

                #output
                (epi_mask,    rest_output, [("brain_mask",          "epi_brain_mask")]),
                (tissue_mask, rest_output, [("out_file",            "tissues_brain_mask")]),
                (realign,     rest_output, [("out_file",            "motion_corrected"),
                                            ("par_file",            "motion_params"),
                                           ]),
                (coreg,       rest_output, [("coregistered_files",  "tissues"),
                                            ("coregistered_source", "anat"),
                                           ]),
                (bandpass,    rest_output, [("out_files",           "time_filtered")]),
                (smooth,      rest_output, [("out_file",            "smooth")]),
              ])

    # add more connections if to perform atlas registration
    if do_atlas:
        coreg_merge  = setup_node(Merge(2), name="coreg_merge")
        wf.connect([
                    (brain_sel,     coreg_merge, [(("out", flatten_list), "in1")]),
                    (rest_input,    coreg_merge, [("atlas_anat",          "in2")]),
                    (coreg_merge,   coreg,       [("out",                 "apply_to_files")]),
                  ])

    else:
        wf.connect([
                    (brain_sel, coreg, [(("out", flatten_list), "apply_to_files")]),
                  ])

    return wf


def attach_rest_preprocessing(main_wf, wf_name="rest_preproc"):
    """ Attach the resting-state MRI pre-processing workflow to the `main_wf`.

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
    # Dependency workflows
    anat_wf  = main_wf.get_node("spm_anat_preproc")
    in_files = get_input_node(main_wf)
    datasink = get_datasink  (main_wf)

    # The workflow box
    rest_wf = rest_preprocessing_wf(wf_name=wf_name)

    # The base name of the 'rest' file for the substitutions
    rest_fbasename = remove_ext(op.basename(get_input_file_name(in_files, 'rest')))

    # prepare substitution for atlas_file, if any
    do_atlas, atlas_file = check_atlas_file()
    if do_atlas:
        atlas_basename = remove_ext(op.basename(atlas_file))

        regexp_subst = [
                         (r"/rw{atlas}\.nii$", "/{atlas}_diff_space.nii"),
                       ]
        regexp_subst = format_pair_list(regexp_subst, atlas=atlas_basename)
        regexp_subst += extension_duplicates(regexp_subst)
        datasink.inputs.regexp_substitutions = extend_trait_list(datasink.inputs.regexp_substitutions,
                                                                 regexp_subst)

    # dataSink output substitutions
    regexp_subst = [
                    (r"/corr_stc{rest}_trim_mean_mask\.\.nii$", "/epi_brain_mask.nii"),
                    (r"/corr_stc{rest}_trim_filt\.nii$",        "/{rest}_time_filt.nii"),
                    (r"/corr_stc{rest}_trim\.nii$",             "/{rest}_motion_corrected.nii"),
                    (r"/stc{rest}_trim\.nii\.par$",             "/motion_parameters.txt"),
                    (r"/rc1[\w]+_corrected\.nii$",              "/coreg_gm.nii"),
                    (r"/rc1[\w]+_corrected_maths\.nii$",        "/tissue_brain_mask.nii"),
                    (r"/rc2[\w]+_corrected\.nii$",              "/coreg_wm.nii"),
                    (r"/rc3[\w]+_corrected\.nii$",              "/coreg_csf.nii"),
                    (r"/rm[\w]+_corrected\.nii$",               "/coreg_anat.nii"),
                   ]
    if regexp_subst:
        regexp_subst  = format_pair_list(regexp_subst, rest=rest_fbasename)
        regexp_subst += extension_duplicates(regexp_subst)
        datasink.inputs.regexp_substitutions = extend_trait_list(datasink.inputs.regexp_substitutions,
                                                                 regexp_subst)

    # input and output anat workflow to main workflow connections
    main_wf.connect([(in_files, rest_wf,  [("rest",                              "rest_input.in_files")]),

                    # anat to fMRI registration
                    (anat_wf,  rest_wf,   [("new_segment.bias_corrected_images", "rest_input.coreg_target"),
                                           ("new_segment.native_class_images",   "rest_input.tissues"),
                                           ("new_segment.bias_corrected_images", "rest_input.anat"),
                                          ]),

                     # test output
                     (rest_wf,  datasink, [
                                           ("rest_output.epi_brain_mask",        "rest.@epi_brain_mask"),
                                           ("rest_output.tissues_brain_mask",    "rest.@tissues_brain_mask"),
                                           ("rest_output.motion_corrected",      "rest.@motion_corr"),
                                           ("rest_output.motion_params",         "rest.@motion_params"),
                                           ("rest_output.tissues",               "rest.@tissues"),
                                           ("rest_output.anat",                  "rest.@anat"),
                                           ("rest_output.time_filtered",         "rest.@time_filtered"),
                                          ]),
                    ])

    if hasattr(anat_wf.inputs, 'atlas_warped'):
            main_wf.connect([(anat_wf,  rest_wf,  [("atlas_warped",           "atlas_anat")]),
                             (rest_wf,  datasink, [("rest_output.atlas",      "rest.@atlas")]),
                            ])

    return main_wf
