# -*- coding: utf-8 -*-
"""
Nipype workflows to process anatomical MRI.
"""
import os.path as op

import nipype.pipeline.engine    as pe
from   nipype.algorithms.misc    import Gunzip
from   nipype.interfaces.utility import Function, Select, IdentityInterface
from   nipype.interfaces         import fsl
from   nipype.interfaces.nipy.preprocess import Trim, ComputeMask

from   .filter   import bandpass_filter
from   .nuisance import rest_noise_filter_wf

from   ..preproc import (auto_spm_slicetime,
                         nipy_motion_correction,
                         spm_coregister,
                         )

from   ..nilearn import mean_img, smooth_img
from   ..nilearn.connectivity import ConnectivityCorrelationInterface
from   ..nilearn.canica import CanICAInterface

from   ..config import setup_node, get_config_setting, check_atlas_file
from   .._utils import format_pair_list, flatten_list
from   ..utils  import (remove_ext,
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

    rest_input.atlas_anat: traits.File
        The atlas file in anatomical space.

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
    rest_output.smooth: traits.File
        The isotropically smoothed time filtered nuisance corrected image.

    rest_output.nuis_corrected: traits.File
        The nuisance corrected fMRI file.

    rest_output.motion_params: traits.File
        The affine transformation file.

    rest_output.time_filtered: traits.File
        The bandpass time filtered fMRI file.

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

    rest_output.motion_regressors: traits.File

    rest_output.compcor_regressors: traits.File

    rest_output.connectivity: traits.File
        If the connectivity node is enabled.
        This is the Numpy matrix text file with the connectivity matrix.

    rest_output.atlas_timeseries: traits.File
        If the connectivity node is enabled.
        This is the Numpy matrix text file with the timeseries extracted using the atlas
        from the resting-state data.

    rest_output.ica_components: traits.File
        If `rest_preproc.canica` is True.

    rest_output.ica_loadings: traits.File
        If `rest_preproc.canica` is True.

    rest_output.ica_score: traits.File
        If `rest_preproc.canica` is True .

    rest_output.all_icc_plot: traits.File
        If `rest_preproc.canica` is True and 'canica_extra.plot' is not False.
        A plot figure in PDF of the ICA results.

    rest_output.iccs_plot: traits.File
        If `rest_preproc.canica` is True and 'canica_extra.plot' is not False.
        A plot figure in PDF of the ICA results.

    Returns
    -------
    wf: nipype Workflow
    """
    # Create the workflow object
    wf = pe.Workflow(name=wf_name)

    # specify input and output fields
    in_fields  = ["in_files",
                  "anat",
                  "atlas_anat",
                  "coreg_target",
                  "tissues",
                  "lowpass_freq",
                  "highpass_freq",]

    out_fields = ["motion_corrected",
                  "motion_params",
                  "tissues",
                  "anat",
                  "time_filtered",
                  "smooth",
                  "tsnr_file",
                  "epi_brain_mask",
                  "tissues_brain_mask",
                  "motion_regressors",
                  "compcor_regressors",
                  "gsr_regressors",
                  "nuis_corrected",]

    do_atlas, _     = check_atlas_file()
    do_connectivity = get_config_setting('rest_preproc.connectivity')
    do_canica       = get_config_setting('rest_preproc.canica')
    do_plot         = get_config_setting('canica_extra.plot', default=True)

    if do_atlas:
        in_fields  += ["atlas_anat"]
        out_fields += ["atlas_rest"]

    if do_atlas and do_connectivity:
        out_fields += ["connectivity", "atlas_timeseries"]

    if do_canica:
        out_fields += ["ica_components", "ica_loadings", "ica_score"]

    if do_canica and do_plot:
        out_fields += ["all_icc_plot", "iccs_plot"]

    # input identities
    rest_input = setup_node(IdentityInterface(fields=in_fields, mandatory_inputs=True),
                            name="rest_input")

    # rs-fMRI preprocessing nodes
    trim    = setup_node(Trim(),
                         name="trim")
    stc_wf  = auto_spm_slicetime()
    realign = setup_node(nipy_motion_correction(), name='realign')

    # average
    average = setup_node(Function(function=mean_img, input_names=["in_file"], output_names=["out_file"],
                                 imports=['from pypes.nilearn import ni2file']),
                        name='average')

    mean_gunzip = setup_node(Gunzip(),        name="mean_gunzip")

    # co-registration nodes
    coreg     = setup_node(spm_coregister(cost_function="mi"), name="coreg_rest")
    brain_sel = setup_node(Select(index=[0, 1, 2]),            name="brain_sel")

    # brain masks
    epi_mask     = setup_node(ComputeMask(),         name='epi_mask')
    tissue_mask  = setup_node(fsl.MultiImageMaths(), name='tissue_mask')
    tissue_mask.inputs.op_string = "-add %s -add %s -kernel gauss 2 -dilM -ero -bin"

    gm_select    = setup_node(Select(index=[0]),     name="gm_sel")
    wmcsf_select = setup_node(Select(index=[1, 2]),  name="wmcsf_sel")

    # noise filter
    noise_wf   = rest_noise_filter_wf()
    wm_select  = setup_node(Select(index=[1]), name="wm_sel")
    csf_select = setup_node(Select(index=[2]), name="csf_sel")

    # bandpass filtering
    bandpass = setup_node(Function(input_names=['files',
                                                'lowpass_freq',
                                                'highpass_freq',
                                                'tr'],
                                   output_names=['out_files'],
                                   function=bandpass_filter),
                          name='bandpass_filter')

    # smooth
    smooth = setup_node(Function(function=smooth_img,
                                 input_names=["in_file", "fwhm"],
                                 output_names=["out_file"],
                                 imports=['from pypes.nilearn import ni2file']),
                         name="fmri_smooth")
    smooth.inputs.fwhm = get_config_setting('fmri_smooth.fwhm', default=8)

    # output identities
    rest_output = setup_node(IdentityInterface(fields=out_fields),
                             name="rest_output")

    # Connect the nodes
    wf.connect([
                # trim
                (rest_input,   trim,         [("in_files", "in_file")]),

                #slice time correction
                (trim,        stc_wf,        [("out_file", "stc_input.in_file")]),

                # motion correction
                (stc_wf,      realign,       [("stc_output.timecorrected_files", "in_file")]),

                # coregistration target
                (realign,     average,       [("out_file", "in_file")]),
                (average,     mean_gunzip,   [("out_file", "in_file")]),
                (mean_gunzip, coreg,         [("out_file", "target")]),

                # epi brain mask
                (mean_gunzip, epi_mask,      [("out_file", "mean_volume")]),

                # coregistration
                (rest_input,  coreg,         [("anat",                "source")]),
                (rest_input,  brain_sel,     [("tissues",             "inlist")]),
                (brain_sel,   coreg,         [(("out", flatten_list), "apply_to_files")]),

                # tissue brain mask
                (coreg,        gm_select,    [("coregistered_files",  "inlist")]),
                (coreg,        wmcsf_select, [("coregistered_files",  "inlist")]),
                (gm_select,    tissue_mask,  [(("out", flatten_list), "in_file")]),
                (wmcsf_select, tissue_mask,  [(("out", flatten_list), "operand_files")]),

                # nuisance correction
                (coreg,         wm_select,  [("coregistered_files",   "inlist",)]),
                (coreg,         csf_select, [("coregistered_files",   "inlist",)]),
                (realign,       noise_wf,   [("out_file",             "rest_noise_input.in_file",)]),
                (tissue_mask,   noise_wf,   [("out_file",             "rest_noise_input.brain_mask")]),
                (wm_select,     noise_wf,   [(("out", flatten_list),  "rest_noise_input.wm_mask")]),
                (csf_select,    noise_wf,   [(("out", flatten_list),  "rest_noise_input.csf_mask")]),
                (realign,       noise_wf,   [("par_file",             "rest_noise_input.motion_params",)]),

                # motion statistics
                # TODO

                # temporal filtering
                (stc_wf,     bandpass,     [("stc_output.time_repetition",          "tr")]),
                (rest_input, bandpass,     [("lowpass_freq",                        "lowpass_freq"),
                                            ("highpass_freq",                       "highpass_freq"),
                                           ]),
                (noise_wf,   bandpass,     [("rest_noise_output.nuis_corrected",    "files")]),

                # smoothing
                (bandpass,    smooth,      [("out_files",           "in_file")]),

                #output
                (epi_mask,    rest_output, [("brain_mask",          "epi_brain_mask")]),
                (tissue_mask, rest_output, [("out_file",            "tissues_brain_mask")]),
                (realign,     rest_output, [("out_file",            "motion_corrected"),
                                            ("par_file",            "motion_params"),
                                           ]),
                (coreg,       rest_output, [("coregistered_files",  "tissues"),
                                            ("coregistered_source", "anat"),
                                           ]),
                (noise_wf,    rest_output, [("rest_noise_output.motion_regressors",  "motion_regressors"),
                                            ("rest_noise_output.compcor_regressors", "compcor_regressors"),
                                            ("rest_noise_output.gsr_regressors",     "gsr_regressors"),
                                            ("rest_noise_output.nuis_corrected",     "nuis_corrected"),
                                            ("rest_noise_output.tsnr_file",          "tsnr_file"),
                                           ]),
                (bandpass,    rest_output, [("out_files",                            "time_filtered")]),
                (smooth,      rest_output, [("out_file",                             "smooth")]),
              ])


    # add more nodes if to perform atlas registration
    if do_atlas:
        coreg_atlas = setup_node(spm_coregister(cost_function="mi"), name="coreg_atlas")

        # set the registration interpolation to nearest neighbour.
        coreg_atlas.inputs.write_interp = 0
        wf.connect([
            (rest_input,  coreg_atlas, [("anat",                "source")]),
            (mean_gunzip, coreg_atlas, [("out_file",            "target")]),
            (rest_input,  coreg_atlas, [("atlas_anat",          "apply_to_files")]),
            (coreg_atlas, rest_output, [("coregistered_files",  "atlas_rest")]),
        ])

    # functional connectivity
    if do_atlas and do_connectivity:
        conn = setup_node(ConnectivityCorrelationInterface(), name="rest_connectivity")
        wf.connect([
            (coreg_atlas, conn,        [("coregistered_files",  "atlas_file")]),
            (bandpass,    conn,        [("out_files",           "in_files")]),
            (conn,        rest_output, [("connectivity",        "connectivity"),
                                        ("timeseries",          "atlas_timeseries"),
                                       ]),
        ])

    # CanICA
    if do_canica:
        ica = setup_node(CanICAInterface(), name="rest_groupica")
        wf.connect([
            (bandpass,    ica,         [("out_files",           "in_files")]),
            (tissue_mask, ica,         [("out_file",            "mask")]),
            (ica,         rest_output, [("components",          "ica_components"),
                                        ("score",               "ica_score"),
                                        ("loadings",            "ica_loadings"),
                                       ]),
        ])

        if do_plot:
            # Connect the plotting nodes
            wf.connect([
                        (ica,         rest_output, [("all_icc_plot", "all_icc_plot"),
                                                    ("iccs_plot",    "iccs_plot"),
                                                   ]),
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

    # dataSink output substitutions
    regexp_subst = [
                    (r"/rc1[\w]+_corrected_maths\.nii$",                          "/tissue_brain_mask.nii"),
                    (r"/rc1[\w]+_corrected\.nii$",                                "/gm.nii"),
                    (r"/rc2[\w]+_corrected\.nii$",                                "/wm.nii"),
                    (r"/rc3[\w]+_corrected\.nii$",                                "/csf.nii"),
                    (r"/rm[\w]+_corrected\.nii$",                                 "/anat.nii"),
                    (r"/corr_stc{rest}_trim\.nii$",                               "/slice_time_corrected.nii"),
                    (r"/stc{rest}_trim\.nii\.par$",                               "/motion_parameters.txt"),
                    (r"/corr_stc{rest}_trim_filt\.nii$",                          "/time_filt.nii"),
                    (r"/corr_stc{rest}_trim_mean_mask\.\.nii$",                   "/epi_brain_mask.nii"),
                    (r"/corr_stc{rest}_trim_filtermotart\.nii$",                  "/motion_corrected.nii"),
                    (r"/corr_stc{rest}_trim_filtermotart[\w_]*_cleaned\.nii$",    "/nuisance_corrected.nii"),
                    (r"/corr_stc{rest}_trim_filtermotart[\w_]*_gsr\.nii$",        "/nuisance_corrected.nii"),
                    (r"/corr_stc{rest}_trim_filtermotart[\w_]*_bandpassed\.nii$", "/time_filtered.nii"),
                    (r"/corr_stc{rest}_trim_filtermotart[\w_]*_smooth\.nii$",     "/smooth.nii"),
                   ]
    regexp_subst = format_pair_list(regexp_subst, rest=rest_fbasename)

    # prepare substitution for atlas_file, if any
    do_atlas, atlas_file = check_atlas_file()
    if do_atlas:
        atlas_basename = remove_ext(op.basename(atlas_file))
        regexp_subst.extend([
                             (r"/[\w]*{atlas}\.nii$", "/{atlas}_{rest}_space.nii"),
                            ])
        regexp_subst = format_pair_list(regexp_subst, atlas=atlas_basename, rest=rest_fbasename)

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
                    (rest_wf,  datasink,  [
                                           ("rest_output.epi_brain_mask",        "rest.@epi_brain_mask"),
                                           ("rest_output.tissues_brain_mask",    "rest.@tissues_brain_mask"),
                                           ("rest_output.tissues",               "rest.@tissues"),
                                           ("rest_output.anat",                  "rest.@anat"),
                                           ("rest_output.motion_regressors",     "rest.@motion_regressors"),
                                           ("rest_output.compcor_regressors",    "rest.@compcor_regressors"),
                                           ("rest_output.gsr_regressors",        "rest.@gsr_regressors"),
                                           ("rest_output.motion_params",         "rest.@motion_params"),
                                           ("rest_output.motion_corrected",      "rest.@motion_corrected"),
                                           ("rest_output.nuis_corrected",        "rest.@nuis_corrected"),
                                           ("rest_output.time_filtered",         "rest.@time_filtered"),
                                           ("rest_output.smooth",                "rest.@smooth"),
                                           ("rest_output.tsnr_file",             "rest.@tsnr"),
                                          ]),
                    ])

    if do_atlas:
            main_wf.connect([(anat_wf,  rest_wf,  [("anat_output.atlas_warped", "rest_input.atlas_anat")]),
                             (rest_wf,  datasink, [("rest_output.atlas_rest",   "rest.@atlas")]),
                            ])

    do_connectivity = get_config_setting('rest_preproc.connectivity')
    if do_atlas and do_connectivity:
            main_wf.connect([(rest_wf,  datasink, [("rest_output.connectivity",     "rest.connectivity.@matrix"),
                                                   ("rest_output.atlas_timeseries", "rest.connectivity.@timeseries"),
                                                  ]),
                            ])

    do_canica = get_config_setting('rest_preproc.canica')
    if do_canica:
            main_wf.connect([(rest_wf,  datasink, [("rest_output.ica_components", "rest.canica.@components"),
                                                   ("rest_output.ica_score",      "rest.canica.@score"),
                                                   ("rest_output.ica_loadings",   "rest.canica.@loadings"),
                                                  ]),
                            ])

    do_plot = get_config_setting('canica_extra.plot', default=True)
    if do_canica and do_plot:
            main_wf.connect([(rest_wf,  datasink, [("rest_output.all_icc_plot", "rest.canica.@all_icc_plot"),
                                                   ("rest_output.iccs_plot",    "rest.canica.@iccs_plot"),
                                                  ]),
                            ])

    return main_wf
