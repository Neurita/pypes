# -*- coding: utf-8 -*-
"""
Nipype workflows to process resting-state functional MRI.
"""
import nipype.pipeline.engine    as pe
import os.path as op
from   nipype.algorithms.misc    import Gunzip
from   nipype.interfaces         import spm, fsl
from   nipype.interfaces.utility import Function, Merge, IdentityInterface

from   .._utils  import format_pair_list, flatten_list
from   ..config  import setup_node
from   ..preproc import (spm_normalize,
                         get_bounding_box,
                         spm_tpm_priors_path,
                         )
from   ..utils   import (remove_ext,
                         selectindex,
                         extend_trait_list,
                         get_input_node,
                         get_datasink,
                         get_input_file_name,
                         extension_duplicates)


def spm_warp_fmri_wf(wf_name="spm_warp_fmri", do_group_template=False):
    """ Run SPM to warp resting-state fMRI pre-processed data to MNI or a given template.

    Tasks:
    - Warping the inputs to MNI or a template, if `do_group_template` is True

    Parameters
    ----------
    wf_name: str

    do_group_template: bool
        If True will expect the wfmri_input.epi_template input and use it as a group template
        for inter-subject registratio.

    Nipype Inputs
    -------------
    wfmri_input.in_file: traits.File
        The slice time and motion corrected fMRI file.

    wfmri_input.time_filtered: traits.File
        The bandpass time filtered fMRI file.

    wfmri_input.avg_epi: traits.File
        The average EPI from the fMRI file.

    wfmri_input.epi_template: traits.File
        Reference EPI template file for inter subject registration.
        If `do_group_template` is True you must specify this input.

    Nipype Outputs
    --------------
    wfmri_output.warped_fmri: traits.File
        The slice time, motion, and nuisance corrected fMRI file registered to the template.

    wfmri_output.wtime_filtered: traits.File
        The bandpass time filtered fMRI file registered to the template.

    wfmri_output.smooth: traits.File
        The smooth bandpass time filtered fMRI file registered to the template.

    wfmri_output.wavg_epi: traits.File
        The average EPI from the fMRI file registered to the template.

    wfmri_output.warp_field: traits.File
        The fMRI to template warp field.

    Returns
    -------
    wf: nipype Workflow
    """
    # Create the workflow object
    wf = pe.Workflow(name=wf_name)

    # specify input and output fields
    in_fields  = ["in_file",
                  "time_filtered",
                  "avg_epi"]

    out_fields = ["warped_fmri",
                  "wtime_filtered",
                  "smooth",
                  "wavg_epi",
                  "warp_field",]

    if do_group_template:
        in_fields += ['epi_template']

    # input identities
    wfmri_input = setup_node(IdentityInterface(fields=in_fields, mandatory_inputs=True),
                             name="wfmri_input")

    # in file unzipper
    in_gunzip = pe.Node(Gunzip(), name="in_gunzip")

    # merge list for normalization input
    merge_list = setup_node(Merge(3), name='merge_for_warp')
    gunzipper = pe.MapNode(Gunzip(), name="gunzip", iterfield=['in_file'])

    # the template bounding box
    tpm_bbox = setup_node(Function(function=get_bounding_box,
                                   input_names=["in_file"],
                                   output_names=["bbox"]),
                          name="tpm_bbox")

    # do template
    if do_group_template:
        gunzip_template = setup_node(Gunzip(), name="gunzip_template",)
        warp = setup_node(spm.Normalize(jobtype="estwrite", out_prefix="wgrptmpl_"),
                          name="fmri_grptemplate_warp",)
        warp_source_arg    = "source"
        warp_outsource_arg = "normalized_source"
        warp_field_arg     = "normalization_parameters"

    else:
        # do MNI
        warp = setup_node(spm_normalize(), name="fmri_warp")
        tpm_bbox.inputs.in_file = spm_tpm_priors_path()
        warp_source_arg    = "image_to_align"
        warp_outsource_arg = "normalized_image"
        warp_field_arg     = "deformation_field"

    # smooth
    # smooth = setup_node(Function(function=smooth_img,
    #                              input_names=["in_file", "fwhm"],
    #                              output_names=["out_file"],
    #                              imports=['from pypes.interfaces.nilearn import ni2file']),
    #                      name="smooth_fmri")
    # smooth.inputs.fwhm = get_config_setting('fmri_smooth.fwhm', default=8)
    # smooth.inputs.out_file = "smooth_{}.nii.gz".format(wf_name)
    smooth = setup_node(fsl.IsotropicSmooth(fwhm=8), name="smooth_fmri")

    # output identities
    rest_output = setup_node(IdentityInterface(fields=out_fields),
                             name="wfmri_output")

    # prepare the inputs
    wf.connect([
                # unzip the in_file input file
                (wfmri_input, in_gunzip,  [("avg_epi", "in_file")]),

                # smooth the time filtered file
                (wfmri_input, smooth,     [("time_filtered", "in_file")]),

                # merge the other input files into a list
                (wfmri_input, merge_list, [("in_file",          "in1"),
                                           ("time_filtered",    "in2"),
                                          ]),
                (smooth,     merge_list, [("out_file",         "in3"), ]),

                # gunzip them for SPM, just in case
                (merge_list, gunzipper, [("out", "in_file")]),

                # bounding box
                (tpm_bbox,  warp,       [("bbox",     "write_bounding_box")]),

                # warp source file
                (in_gunzip, warp,       [("out_file",  warp_source_arg)]),

                # apply to files
                (gunzipper, warp,       [("out_file", "apply_to_files")]),

                # outputs
                (warp, rest_output,     [(warp_field_arg,     "warp_field"),
                                         (warp_outsource_arg, "wavg_epi"),
                                        ]),

               ])

    if do_group_template:
        wf.connect([
                    # unzip and forward the template file
                    (wfmri_input,      gunzip_template, [("epi_template", "in_file")]),
                    (gunzip_template, warp,             [("out_file",     "template")]),

                    # get template bounding box to apply to results
                    (wfmri_input,      tpm_bbox,        [("epi_template", "in_file")]),
        ])

    wf.connect([
                # smooth
                #(warp,   smooth,      [(("normalized_files", selectindex, [1]), "in_file")]),

                # output
                (warp,   rest_output, [(("normalized_files", selectindex, [0]), "warped_fmri"),
                                       (("normalized_files", selectindex, [1]), "wtime_filtered"),
                                       (("normalized_files", selectindex, [2]), "smooth"),
                                      ]),
               ])

    return wf


def attach_spm_warp_fmri_wf(main_wf, registration_wf_name="spm_warp_fmri", do_group_template=False):
    """ Attach the fMRI inter-subject spatial normalization workflow to the `main_wf`.

    Parameters
    ----------
    main_wf: nipype Workflow

    registration_wf_name: str
        Name of the registration workflow.

    do_group_template: bool
        If True will attach the group template creation and pre-processing pipeline.

    Nipype Inputs for `main_wf`
    ---------------------------
    Note: The `main_wf` workflow is expected to have an `input_files` and a `datasink` nodes.

    input_files.select.anat: input node

    datasink: nipype Node

    Workflow Dependencies
    ---------------------
    fmri_cleanup

    spm_fmri_template, if do_group_template is True

    Returns
    -------
    main_wf: nipype Workflow
    """
    # Dependency workflows
    cleanup_wf = main_wf.get_node("fmri_cleanup")

    in_files = get_input_node(main_wf)
    datasink = get_datasink  (main_wf)

    if do_group_template:
        template_name = 'grptemplate'
    else:
        template_name = 'mni'

    warp_wf_name = "{}_{}".format(registration_wf_name, template_name)
    warp_fmri_wf = spm_warp_fmri_wf(warp_wf_name, do_group_template=do_group_template)

    # dataSink output substitutions
    # The base name of the 'rest' file for the substitutions
    rest_fbasename = remove_ext(op.basename(get_input_file_name(in_files, 'rest')))
    anat_fbasename = remove_ext(op.basename(get_input_file_name(in_files, 'anat')))

    regexp_subst = [
                    (r"/corr_stc{fmri}_trim_mean_sn.mat$", "/{fmri}_grptemplate_params.mat"),
                    (r"/y_corr_stc{fmri}_trim_mean\.nii$", "/{fmri}_to_mni_warpfield.nii"),

                    (r"/wgrptmpl_corr_stc{fmri}_trim_mean\.nii$",                          "/avg_epi_grptemplate.nii"),
                    (r"/wgrptmpl_corr_stc{fmri}_trim\.nii$",                               "/{fmri}_trimmed_grptemplate.nii"),
                    (r"/wgrptmpl_corr_stc{fmri}_trim_filtermotart[\w_]*_cleaned\.nii$",    "/{fmri}_nuisance_corrected_grptemplate.nii"),
                    (r"/wgrptmpl_corr_stc{fmri}_trim_filtermotart[\w_]*_gsr\.nii$",        "/{fmri}_nuisance_corrected_grptemplate.nii"),
                    (r"/wgrptmpl_corr_stc{fmri}_trim_filtermotart[\w_]*_bandpassed\.nii$", "/{fmri}_time_filtered_grptemplate.nii"),
                    (r"/wgrptmpl_corr_stc{fmri}_trim_filtermotart[\w_]*_smooth\.nii$",     "/{fmri}_smooth_grptemplate.nii"),

                    (r"/wcorr_stc{fmri}_trim_mean\.nii$",                          "/avg_epi_mni.nii"),
                    (r"/wcorr_stc{fmri}_trim\.nii$",                               "/{fmri}_trimmed_mni.nii"),
                    (r"/wcorr_stc{fmri}_trim_filtermotart[\w_]*_cleaned\.nii$",    "/{fmri}_nuisance_corrected_mni.nii"),
                    (r"/wcorr_stc{fmri}_trim_filtermotart[\w_]*_gsr\.nii$",        "/{fmri}_nuisance_corrected_mni.nii"),
                    (r"/wcorr_stc{fmri}_trim_filtermotart[\w_]*_bandpassed\.nii$", "/{fmri}_time_filtered_mni.nii"),
                    (r"/wcorr_stc{fmri}_trim_filtermotart[\w_]*_smooth\.nii$",     "/{fmri}_smooth_mni.nii"),
                   ]
    regexp_subst = format_pair_list(regexp_subst, fmri=rest_fbasename, anat=anat_fbasename)
    regexp_subst += extension_duplicates(regexp_subst)
    datasink.inputs.regexp_substitutions = extend_trait_list(datasink.inputs.regexp_substitutions,
                                                             regexp_subst)

    # input and output anat workflow to main workflow connections
    main_wf.connect([# clean_up_wf to registration_wf
                     (cleanup_wf, warp_fmri_wf, [
                                                 ("rest_output.motion_corrected", "wfmri_input.in_file"),
                                                 ("rest_output.time_filtered",    "wfmri_input.time_filtered"),
                                                 ("rest_output.avg_epi",          "wfmri_input.avg_epi"),
                                                ]),

                    # output
                    (warp_fmri_wf,  datasink,  [
                                                ("wfmri_output.warped_fmri",     "rest.@warped_fmri_{}".format(template_name)),
                                                ("wfmri_output.wtime_filtered",  "rest.@time_filtered_{}".format(template_name)),
                                                ("wfmri_output.smooth",          "rest.@smooth_{}".format(template_name)),
                                                ("wfmri_output.wavg_epi",        "rest.@avg_epi_{}".format(template_name)),
                                                ("wfmri_output.warp_field",      "rest.@warp_field_{}".format(template_name)),
                                               ]),
                    ])

    return main_wf
