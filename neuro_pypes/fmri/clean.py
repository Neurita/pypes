# -*- coding: utf-8 -*-
"""
Nipype workflow to clean up resting-state functional MRI.
"""
import nipype.pipeline.engine as pe
import os.path as op
from   nipype.algorithms.misc import Gunzip
from   nipype.interfaces import fsl
from   nipype.interfaces.nipy.preprocess import Trim, ComputeMask
from   nipype.interfaces.utility import Function, Select, IdentityInterface
from   neuro_pypes.interfaces.nilearn import mean_img, smooth_img

from   .filter import bandpass_filter
from   .nuisance import rest_noise_filter_wf
from   .._utils import format_pair_list, flatten_list
from   ..config import setup_node, get_config_setting
from   ..preproc import (auto_spm_slicetime,
                         nipy_motion_correction,
                         spm_coregister,
                         )
from   ..utils  import (remove_ext,
                        extend_trait_list,
                        get_input_node,
                        get_interface_node,
                        get_datasink,
                        get_input_file_name,
                        extension_duplicates)


def fmri_cleanup_wf(wf_name="fmri_cleanup"):
    """ Run the resting-state fMRI pre-processing workflow against the rest files in `data_dir`.

    Tasks:
    - Trim first 6 volumes of the rs-fMRI file.
    - Slice Timing correction.
    - Motion and nuisance correction.
    - Calculate brain mask in fMRI space.
    - Bandpass frequency filtering for resting-state fMRI.
    - Smoothing.
    - Tissue maps co-registration to fMRI space.

    Parameters
    ----------
    wf_name: str

    Nipype Inputs
    -------------
    rest_input.in_file: traits.File
        The resting-state fMRI file.

    rest_input.anat: traits.File
        Path to the high-contrast anatomical image.

    rest_input.tissues: list of traits.File
        Paths to the tissue segmentations in anatomical space.
        Expected to have this order: GM, WM and CSF.

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

    rest_output.tissues: list of traits.File
        The tissues segmentation volume in fMRI space.
        Expected to have this order: GM, WM and CSF.

    rest_output.anat: traits.File
        The T1w image in fMRI space.

    rest_output.avg_epi: traits.File
        The average EPI image in fMRI space after slice-time and motion correction.

    rest_output.motion_regressors: traits.File

    rest_output.compcor_regressors: traits.File

    rest_output.art_displacement_files
        One image file containing the voxel-displacement timeseries.

    rest_output.art_intensity_files
        One file containing the global intensity values determined from the brainmask.

    rest_output.art_norm_files
        One file containing the composite norm.

    rest_output.art_outlier_files
         One file containing a list of 0-based indices corresponding to outlier volumes.

    rest_output.art_plot_files
        One image file containing the detected outliers.

    rest_output.art_statistic_files
        One file containing information about the different types of artifacts and if design info is provided then
        details of stimulus correlated motion and a listing or artifacts by event type.

    Returns
    -------
    wf: nipype Workflow
    """
    # Create the workflow object
    wf = pe.Workflow(name=wf_name)

    # specify input and output fields
    in_fields  = ["in_file",
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
                  "avg_epi",
                  "time_filtered",
                  "smooth",
                  "tsnr_file",
                  "epi_brain_mask",
                  "tissues_brain_mask",
                  "motion_regressors",
                  "compcor_regressors",
                  "gsr_regressors",
                  "nuis_corrected",
                  "art_displacement_files",
                  "art_intensity_files",
                  "art_norm_files",
                  "art_outlier_files",
                  "art_plot_files",
                  "art_statistic_files",]

    # input identities
    rest_input = setup_node(IdentityInterface(fields=in_fields, mandatory_inputs=True),
                            name="rest_input")

    # rs-fMRI preprocessing nodes
    trim    = setup_node(Trim(), name="trim")

    stc_wf  = auto_spm_slicetime()
    realign = setup_node(nipy_motion_correction(), name='realign')

    # average
    average = setup_node(Function(function=mean_img, input_names=["in_file"], output_names=["out_file"],
                                  imports=['from neuro_pypes.interfaces.nilearn import ni2file']),
                         name='average_epi')

    mean_gunzip = setup_node(Gunzip(), name="mean_gunzip")

    # co-registration nodes
    coreg     = setup_node(spm_coregister(cost_function="mi"), name="coreg_fmri")
    brain_sel = setup_node(Select(index=[0, 1, 2]),            name="brain_sel")

    # brain mask made with EPI
    epi_mask    = setup_node(ComputeMask(),         name='epi_mask')

    # brain mask made with the merge of the tissue segmentations
    tissue_mask = setup_node(fsl.MultiImageMaths(), name='tissue_mask')
    tissue_mask.inputs.op_string = "-add %s -add %s -abs -kernel gauss 4 -dilM -ero -kernel gauss 1 -dilM -bin"
    tissue_mask.inputs.out_file = "tissue_brain_mask.nii.gz"

    # select tissues
    gm_select    = setup_node(Select(index=[0]),    name="gm_sel")
    wmcsf_select = setup_node(Select(index=[1, 2]), name="wmcsf_sel")

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
                          name='bandpass')

    # smooth
    smooth = setup_node(Function(function=smooth_img,
                                 input_names=["in_file", "fwhm"],
                                 output_names=["out_file"],
                                 imports=['from neuro_pypes.interfaces.nilearn import ni2file']),
                         name="smooth")
    smooth.inputs.fwhm = get_config_setting('fmri_smooth.fwhm', default=8)
    smooth.inputs.out_file = "smooth_std_{}.nii.gz".format(wf_name)

    # output identities
    rest_output = setup_node(IdentityInterface(fields=out_fields),
                             name="rest_output")

    # Connect the nodes
    wf.connect([
                # trim
                (rest_input,   trim,         [("in_file", "in_file")]),

                #slice time correction
                (trim,        stc_wf,        [("out_file", "stc_input.in_file")]),

                # motion correction
                (stc_wf,      realign,       [("stc_output.timecorrected_files", "in_file")]),

                # coregistration target
                (realign,     average,       [("out_file", "in_file")]),
                (average,     mean_gunzip,   [("out_file", "in_file")]),
                (mean_gunzip, coreg,         [("out_file", "target")]),

                # epi brain mask
                (average,     epi_mask,      [("out_file", "mean_volume")]),

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

                # temporal filtering
                (noise_wf,    bandpass,    [("rest_noise_output.nuis_corrected", "files")]),
                #(realign,     bandpass,    [("out_file", "files")]),
                (stc_wf,      bandpass,    [("stc_output.time_repetition", "tr")]),
                (rest_input,  bandpass,    [("lowpass_freq",               "lowpass_freq"),
                                            ("highpass_freq",              "highpass_freq"),
                                           ]),
                (bandpass,     smooth,     [("out_files", "in_file")]),

                # output
                (epi_mask,    rest_output, [("brain_mask", "epi_brain_mask")]),
                (tissue_mask, rest_output, [("out_file",   "tissues_brain_mask")]),
                (realign,     rest_output, [("out_file",   "motion_corrected"),
                                            ("par_file",   "motion_params"),
                                           ]),
                (coreg,       rest_output, [("coregistered_files",  "tissues"),
                                            ("coregistered_source", "anat"),
                                           ]),
                (noise_wf,    rest_output, [("rest_noise_output.motion_regressors",      "motion_regressors"),
                                            ("rest_noise_output.compcor_regressors",     "compcor_regressors"),
                                            ("rest_noise_output.gsr_regressors",         "gsr_regressors"),
                                            ("rest_noise_output.nuis_corrected",         "nuis_corrected"),
                                            ("rest_noise_output.tsnr_file",              "tsnr_file"),
                                            ("rest_noise_output.art_displacement_files", "art_displacement_files"),
                                            ("rest_noise_output.art_intensity_files",    "art_intensity_files"),
                                            ("rest_noise_output.art_norm_files",         "art_norm_files"),
                                            ("rest_noise_output.art_outlier_files",      "art_outlier_files"),
                                            ("rest_noise_output.art_plot_files",         "art_plot_files"),
                                            ("rest_noise_output.art_statistic_files",    "art_statistic_files"),
                                           ]),
                (average,     rest_output, [("out_file",  "avg_epi")]),
                (bandpass,    rest_output, [("out_files", "time_filtered")]),
                (smooth,      rest_output, [("out_file",  "smooth")]),
              ])

    return wf


def attach_fmri_cleanup_wf(main_wf, wf_name="fmri_cleanup"):
    """ Attach the resting-state MRI pre-processing workflow to the `main_wf`.

    Parameters
    ----------
    main_wf: nipype Workflow

    registration_wf_name: str
        Name of the registration workflow.

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
    in_files = get_input_node(main_wf)
    datasink = get_datasink  (main_wf)

    anat_output = get_interface_node(main_wf, "anat_output")

    # create the fMRI preprocessing pipelines
    cleanup_wf = fmri_cleanup_wf(wf_name)

    # dataSink output substitutions
    # The base name of the 'rest' file for the substitutions
    rest_fbasename = remove_ext(op.basename(get_input_file_name(in_files, 'rest')))
    anat_fbasename = remove_ext(op.basename(get_input_file_name(in_files, 'anat')))

    regexp_subst = [
                    (r"/rc1[\w]+_corrected\.nii$",                                 "/gm_{rest}.nii"),
                    (r"/rc2[\w]+_corrected\.nii$",                                 "/wm_{rest}.nii"),
                    (r"/rc3[\w]+_corrected\.nii$",                                 "/csf_{rest}.nii"),
                    (r"/rm[\w]+_corrected\.nii$",                                  "/{anat}_{rest}.nii"),
                    (r"/corr_stc{rest}_trim\.nii$",                                "/slice_time_corrected.nii"),
                    (r"/stc{rest}_trim\.nii\.par$",                                "/motion_parameters.txt"),
                    (r"/corr_stc{rest}_trim_filt\.nii$",                           "/time_filt.nii"),
                    (r"/corr_stc{rest}_trim_mean_mask\.\.nii$",                    "/epi_brain_mask_{rest}.nii"),
                    (r"/tissue_brain_mask\.nii$",                                  "/tissue_brain_mask_{rest}.nii"),
                    (r"/corr_stc{rest}_trim_mean\.nii$",                           "/avg_epi.nii"),

                    (r"/art\..*_outliers\.txt$",                                   "/artifact_outliers.txt"),
                    (r"/global_intensity\..*\.txt$",                               "/global_intensities.txt"),
                    (r"/norm\..*_outliers\.txt$",                                  "/motion_norms.txt"),
                    (r"/stats\..*\.txt$",                                          "/motion_stats.json"),
                    (r"/plot\..*\.png$",                                           "/artifact_plots.png"),

                    (r"/corr_stc{rest}_trim_filtermotart\.nii$",                   "/{rest}_motion_corrected.nii"),
                    (r"/corr_stc{rest}_trim_filtermotart[\w_]*_cleaned\.nii$",     "/{rest}_nuisance_corrected.nii"),
                    (r"/corr_stc{rest}_trim_filtermotart[\w_]*_gsr\.nii$",         "/{rest}_nuisance_corrected.nii"),
                    (r"/corr_stc{rest}_trim_filtermotart[\w_]*_bandpassed\.nii$",  "/{rest}_time_filtered.nii"),
                    (r"/corr_stc{rest}_trim_filtermotart[\w_]*_smooth\.nii$",      "/{rest}_smooth.nii"),
                   ]
    regexp_subst = format_pair_list(regexp_subst, rest=rest_fbasename, anat=anat_fbasename)
    regexp_subst += extension_duplicates(regexp_subst)
    datasink.inputs.regexp_substitutions = extend_trait_list(datasink.inputs.regexp_substitutions,
                                                             regexp_subst)

    # input and output anat workflow to main workflow connections
    main_wf.connect([(in_files,   cleanup_wf, [("rest", "rest_input.in_file")]),

                    # anat to fMRI registration inputs
                    (anat_output, cleanup_wf, [("tissues_native", "rest_input.tissues"),
                                               ("anat_biascorr",  "rest_input.anat"),
                                              ]),

                    # clean_up_wf to datasink
                    (cleanup_wf,  datasink,  [
                                                ("rest_output.epi_brain_mask",         "rest.@epi_brain_mask"),
                                                ("rest_output.tissues_brain_mask",     "rest.@tissues_brain_mask"),
                                                ("rest_output.tissues",                "rest.@tissues"),
                                                ("rest_output.anat",                   "rest.@anat"),
                                                ("rest_output.motion_regressors",      "rest.@motion_regressors"),
                                                ("rest_output.compcor_regressors",     "rest.@compcor_regressors"),
                                                ("rest_output.gsr_regressors",         "rest.@gsr_regressors"),
                                                ("rest_output.motion_params",          "rest.@motion_params"),
                                                ("rest_output.motion_corrected",       "rest.@motion_corrected"),
                                                ("rest_output.nuis_corrected",         "rest.@nuis_corrected"),
                                                ("rest_output.time_filtered",          "rest.@time_filtered"),
                                                ("rest_output.smooth",                 "rest.@smooth"),
                                                ("rest_output.avg_epi",                "rest.@avg_epi"),
                                                ("rest_output.tsnr_file",              "rest.@tsnr"),
                                                ("rest_output.art_displacement_files", "rest.artifact_stats.@displacement"),
                                                ("rest_output.art_intensity_files",    "rest.artifact_stats.@art_intensity"),
                                                ("rest_output.art_norm_files",         "rest.artifact_stats.@art_norm"),
                                                ("rest_output.art_outlier_files",      "rest.artifact_stats.@art_outlier"),
                                                ("rest_output.art_plot_files",         "rest.artifact_stats.@art_plot"),
                                                ("rest_output.art_statistic_files",    "rest.artifact_stats.@art_statistic"),
                                               ]),
                    ])

    return main_wf
