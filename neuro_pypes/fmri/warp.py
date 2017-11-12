# -*- coding: utf-8 -*-
"""
Nipype workflows to process resting-state functional MRI.
"""
import os.path as op

import nipype.pipeline.engine    as pe
from   nipype.algorithms.misc    import Gunzip
from   nipype.interfaces.utility import Function, Merge, IdentityInterface
from   nipype.interfaces import spm, fsl

from   .._utils  import format_pair_list
from   ..config  import setup_node, get_config_setting, check_atlas_file
from   ..preproc import (spm_normalize,
                         get_bounding_box,
                         spm_tpm_priors_path,
                         spm_coregister,
                         spm_apply_deformations,
                         )
from   ..utils   import (remove_ext,
                         selectindex,
                         extend_trait_list,
                         get_input_node,
                         get_datasink,
                         get_input_file_name,
                         extension_duplicates)


def spm_warp_fmri_wf(wf_name="spm_warp_fmri", register_to_grptemplate=False):
    """ Run SPM to warp resting-state fMRI pre-processed data to MNI or a given 
    template.

    Tasks:
    - Warping the inputs to MNI or a template, if `do_group_template` is True

    Parameters
    ----------
    wf_name: str

    register_to_grptemplate: bool
        If True will expect the wfmri_input.epi_template input and use it as a group template
        for inter-subject registratio.

    Nipype Inputs
    -------------
    wfmri_input.in_file: traits.File
        The slice time and motion corrected fMRI file.

    wfmri_input.reference_file: traits.File
        The anatomical image in its native space
        for registration reference.

    wfmri_input.anat_fmri: traits.File
        The anatomical image in fMRI space.

    wfmri_input.anat_to_mni_warp: traits.File
        The warp field from the transformation of the
        anatomical image to the standard MNI space.

    wfmri_input.time_filtered: traits.File
        The bandpass time filtered fMRI file.

    wfmri_input.avg_epi: traits.File
        The average EPI from the fMRI file.

    wfmri_input.epi_template: traits.File
        Reference EPI template file for inter subject registration.
        If `do_group_template` is True you must specify this input.

    wfmri_input.brain_mask: traits.File
        Brain mask in fMRI space.

    wfmri_input.atlas_anat: traits.File
        Atlas in subject anatomical space.

    Nipype Outputs
    --------------
    wfmri_output.warped_fmri: traits.File
        The slice time, motion, and nuisance corrected fMRI
        file registered to the template.

    wfmri_output.wtime_filtered: traits.File
        The bandpass time filtered fMRI file
        registered to the template.

    wfmri_output.smooth: traits.File
        The smooth bandpass time filtered fMRI file
        registered to the template.

    wfmri_output.wavg_epi: traits.File
        The average EPI from the fMRI file
        registered to the template.

    wfmri_output.warp_field: traits.File
        The fMRI to template warp field.

    wfmri_output.coreg_avg_epi: traits.File
        The average EPI image in anatomical space.

        Only if registration.anat2fmri is false.

    wfmri_output.coreg_others: traits.File
        Other mid-preprocessing fmri images registered to
        anatomical space:

        - wfmri_input.in_file,

        - wfmri_input.brain_mask,

        - wfmri_input.time_filtered.

        Only if registration.anat2fmri is false

    wfmri_output.wbrain_mask: traits.File
        Brain mask in fMRI space warped to MNI.

    Returns
    -------
    wf: nipype Workflow
    """
    # Create the workflow object
    wf = pe.Workflow(name=wf_name)

    # specify input and output fields
    in_fields  = ["in_file",
                  "anat_fmri",
                  "anat_to_mni_warp",
                  "brain_mask",
                  "reference_file",
                  "time_filtered",
                  "avg_epi",]

    out_fields = ["warped_fmri",
                  "wtime_filtered",
                  "smooth",
                  "wavg_epi",
                  "wbrain_mask",
                  "warp_field",
                  "coreg_avg_epi",
                  "coreg_others"
                  ]

    if register_to_grptemplate:
        in_fields += ['epi_template']

    do_atlas, _ = check_atlas_file()
    if do_atlas:
        in_fields  += ["atlas_anat"]
        out_fields += ["atlas_fmri"]

    # input identities
    wfmri_input = setup_node(IdentityInterface(fields=in_fields, mandatory_inputs=True),
                             name="wfmri_input")

    # in file unzipper
    in_gunzip = pe.Node(Gunzip(), name="in_gunzip")

    # merge list for normalization input
    merge_list = pe.Node(Merge(2), name='merge_for_warp')
    gunzipper = pe.MapNode(Gunzip(), name="gunzip", iterfield=['in_file'])

    # the template bounding box
    tpm_bbox = setup_node(Function(function=get_bounding_box,
                                   input_names=["in_file"],
                                   output_names=["bbox"]),
                          name="tpm_bbox")

    # smooth the final result
    smooth = setup_node(fsl.IsotropicSmooth(fwhm=8, output_type='NIFTI'), name="smooth_fmri")

    # output identities
    rest_output = setup_node(IdentityInterface(fields=out_fields),
                             name="wfmri_output")


    # check how to perform the registration, to decide how to build the pipeline
    anat2fmri = get_config_setting('registration.anat2fmri', False)
    # register to group template
    if register_to_grptemplate:
        gunzip_template = pe.Node(Gunzip(), name="gunzip_template",)
        warp = setup_node(spm.Normalize(jobtype="estwrite", out_prefix="wgrptmpl_"),
                          name="fmri_grptemplate_warp",)
        warp_source_arg    = "source"
        warp_outsource_arg = "normalized_source"
        warp_field_arg     = "normalization_parameters"

    elif anat2fmri:
        # register to standard template
        warp = setup_node(spm_normalize(), name="fmri_warp")
        tpm_bbox.inputs.in_file = spm_tpm_priors_path()
        warp_source_arg    = "image_to_align"
        warp_outsource_arg = "normalized_image"
        warp_field_arg     = "deformation_field"

    else: # anat2fmri is False
        coreg       = setup_node(spm_coregister(cost_function="mi"), name="coreg_fmri")
        warp        = setup_node(spm_apply_deformations(), name="fmri_warp")
        coreg_files = pe.Node(Merge(3), name='merge_for_coreg')
        warp_files  = pe.Node(Merge(2), name='merge_for_warp')
        tpm_bbox.inputs.in_file = spm_tpm_priors_path()

    # make the connections
    if register_to_grptemplate:
        wf.connect([
                    # get template bounding box to apply to results
                    (wfmri_input, tpm_bbox, [("epi_template", "in_file")]),

                    # unzip and forward the template file
                    (wfmri_input,     gunzip_template, [("epi_template", "in_file")]),
                    (gunzip_template, warp,            [("out_file",     "template")]),

                    # get template bounding box to apply to results
                    (wfmri_input, tpm_bbox, [("epi_template", "in_file")]),
                    ])

    if anat2fmri or register_to_grptemplate:
        # prepare the inputs
        wf.connect([
                    # unzip the in_file input file
                    (wfmri_input, in_gunzip, [("avg_epi", "in_file")]),

                    # warp source file
                    (in_gunzip, warp, [("out_file", warp_source_arg)]),

                    # bounding box
                    (tpm_bbox,  warp, [("bbox", "write_bounding_box")]),

                    # merge the other input files into a list
                    (wfmri_input, merge_list, [("in_file",       "in1"),
                                               ("time_filtered", "in2"),
                                              ]),

                    # gunzip them for SPM
                    (merge_list, gunzipper, [("out", "in_file")]),

                    # apply to files
                    (gunzipper, warp,       [("out_file", "apply_to_files")]),

                    # outputs
                    (warp, rest_output,     [(warp_field_arg,     "warp_field"),
                                             (warp_outsource_arg, "wavg_epi"),
                                            ]),

                   ])

    else: # FMRI to ANAT
        wf.connect([
                    (wfmri_input, coreg,      [("reference_file", "target")]),

                    # unzip the in_file input file
                    (wfmri_input, in_gunzip,  [("avg_epi",  "in_file")]),
                    (in_gunzip,   coreg,      [("out_file", "source")]),

                    # merge the other input files into a list
                    (wfmri_input, coreg_files, [("in_file",       "in1"),
                                                ("time_filtered", "in2"),
                                                ("brain_mask",    "in3"),
                                               ]),

                    # gunzip them for SPM
                    (coreg_files, gunzipper, [("out", "in_file")]),

                    # coregister fmri to anat
                    (gunzipper,   coreg,     [("out_file", "apply_to_files")]),

                    # anat to mni warp field
                    (wfmri_input, warp, [("anat_to_mni_warp", "deformation_file")]),

                    # bounding box
                    (tpm_bbox,  warp, [("bbox", "write_bounding_box")]),

                    # apply to files
                    (coreg, warp_files, [("coregistered_source", "in1")]),
                    (coreg, warp_files, [("coregistered_files",  "in2")]),

                    (warp_files, warp, [("out", "apply_to_files")]),

                    # outputs
                    (warp, rest_output,  [("normalized_files",  "warped_files"),]),
                    (warp, rest_output,  [(("normalized_files", selectindex, 0), "wavg_epi"),]),

                    (coreg, rest_output, [("coregistered_source", "coreg_avg_epi")]),
                    (coreg, rest_output, [("coregistered_files",  "coreg_others")]),
                   ])

    # atlas file in fMRI space
    if anat2fmri:
        coreg_atlas = setup_node(spm_coregister(cost_function="mi"), name="coreg_atlas2fmri")

        # set the registration interpolation to nearest neighbour.
        coreg_atlas.inputs.write_interp = 0
        wf.connect([
                    (wfmri_input, coreg_atlas, [("reference_file", "source"),
                                                ("atlas_anat", "apply_to_files"),
                                               ]),
                    (in_gunzip,   coreg_atlas, [("out_file",           "target")]),
                    (coreg_atlas, rest_output, [("coregistered_files", "atlas_fmri")]),
                  ])

    # smooth and sink
    wf.connect([
                # smooth the final bandpassed image
                (warp,   smooth,      [(("normalized_files", selectindex, 1), "in_file")]),

                # output
                (smooth, rest_output, [("out_file", "smooth")]),
                (warp,   rest_output, [(("normalized_files", selectindex, 0), "warped_fmri"),
                                       (("normalized_files", selectindex, 1), "wtime_filtered"),
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
    fmri_cleanup, the cleanup and preprocessing of the fMRI data

    spm_anat_preproc, for the anatomical to MNI space transformation

    spm_fmri_template, if do_group_template is True

    Returns
    -------
    main_wf: nipype Workflow
    """
    # Dependency workflows
    anat_wf    = main_wf.get_node("spm_anat_preproc")
    cleanup_wf = main_wf.get_node("fmri_cleanup")

    in_files = get_input_node(main_wf)
    datasink = get_datasink  (main_wf)

    if do_group_template:
        template_name = 'grptemplate'
    else:
        template_name = 'stdtemplate'

    warp_wf_name = "{}_{}".format(registration_wf_name, template_name)
    warp_fmri_wf = spm_warp_fmri_wf(warp_wf_name, register_to_grptemplate=do_group_template)

    # dataSink output substitutions
    # The base name of the 'rest' file for the substitutions
    rest_fbasename = remove_ext(op.basename(get_input_file_name(in_files, 'rest')))
    anat_fbasename = remove_ext(op.basename(get_input_file_name(in_files, 'anat')))

    regexp_subst = [
                    (r"/corr_stc{fmri}_trim_mean_sn.mat$", "/{fmri}_grptemplate_params.mat"),
                    (r"/y_corr_stc{fmri}_trim_mean\.nii$", "/{fmri}_to_mni_warpfield.nii"),
                    (r"/rcorr_stc{fmri}_trim_mean.nii$",   "/avg_epi_anat.nii"),

                    (r"/wgrptmpl_corr_stc{fmri}_trim_mean\.nii$",                          "/avg_epi_grptemplate.nii"),
                    (r"/wgrptmpl_corr_stc{fmri}_trim\.nii$",                               "/{fmri}_trimmed_grptemplate.nii"),
                    (r"/wgrptmpl_corr_stc{fmri}_trim_filtermotart[\w_]*_cleaned\.nii$",    "/{fmri}_nuisance_corrected_grptemplate.nii"),
                    (r"/wgrptmpl_corr_stc{fmri}_trim_filtermotart[\w_]*_gsr\.nii$",        "/{fmri}_nuisance_corrected_grptemplate.nii"),
                    (r"/wgrptmpl_corr_stc{fmri}_trim_filtermotart[\w_]*_bandpassed\.nii$", "/{fmri}_time_filtered_grptemplate.nii"),
                    (r"/wgrptmpl_corr_stc{fmri}_trim_filtermotart[\w_]*_smooth\.nii$",     "/{fmri}_smooth_grptemplate.nii"),

                    (r"/w[r]?corr_stc{fmri}_trim_mean\.nii$",                          "/avg_epi_mni.nii"),
                    (r"/w[r]?corr_stc{fmri}_trim\.nii$",                               "/{fmri}_trimmed_mni.nii"),
                    (r"/w[r]?corr_stc{fmri}_trim_filtermotart[\w_]*_cleaned\.nii$",    "/{fmri}_nuisance_corrected_mni.nii"),
                    (r"/w[r]?corr_stc{fmri}_trim_filtermotart[\w_]*_gsr\.nii$",        "/{fmri}_nuisance_corrected_mni.nii"),
                    (r"/w[r]?corr_stc{fmri}_trim_filtermotart[\w_]*_bandpassed\.nii$", "/{fmri}_time_filtered_mni.nii"),
                    (r"/w[r]?corr_stc{fmri}_trim_filtermotart[\w_]*_smooth\.nii$",     "/{fmri}_smooth_mni.nii"),

                    (r"/w[r]?corr_stc{fmri}_trim[\w_]*_smooth\.nii$",       "/{fmri}_nofilt_smooth_mni.nii"),
                    (r"/w[r]?corr_stc{fmri}_trim[\w_]*_cleaned\.nii$",      "/{fmri}_nofilt_nuisance_corrected_mni.nii"),
                    (r"/w[r]?corr_stc{fmri}_trim[\w_]*_gsr\.nii$",          "/{fmri}_nofilt_nuisance_corrected_mni.nii"),
                    (r"/w[r]?corr_stc{fmri}_trim[\w_]*_bandpassed\.nii$",   "/{fmri}_nofilt_time_filtered_mni.nii"),
                    (r"/w[r]?corr_stc{fmri}_trim[\w_]*_smooth\.nii$",       "/{fmri}_nofilt_smooth_mni.nii"),
                  ]
    regexp_subst = format_pair_list(regexp_subst, fmri=rest_fbasename, anat=anat_fbasename)

    # prepare substitution for atlas_file, if any
    do_atlas, atlas_file = check_atlas_file()
    if do_atlas:
        atlas_basename = remove_ext(op.basename(atlas_file))
        regexp_subst.extend([
                             (r"/[\w]*{atlas}.*\.nii$", "/{atlas}_{fmri}_space.nii"),
                            ])
        regexp_subst = format_pair_list(regexp_subst, atlas=atlas_basename,
                                                      fmri=rest_fbasename)

    regexp_subst += extension_duplicates(regexp_subst)
    datasink.inputs.regexp_substitutions = extend_trait_list(datasink.inputs.regexp_substitutions,
                                                             regexp_subst)

    # input and output anat workflow to main workflow connections
    main_wf.connect([# clean_up_wf to registration_wf
                     (cleanup_wf, warp_fmri_wf, [
                                                 ("rest_output.motion_corrected",   "wfmri_input.in_file"),
                                                 ("rest_output.anat",               "wfmri_input.anat_fmri"),
                                                 ("rest_output.time_filtered",      "wfmri_input.time_filtered"),
                                                 ("rest_output.avg_epi",            "wfmri_input.avg_epi"),
                                                 ("rest_output.tissues_brain_mask", "wfmri_input.brain_mask"),
                                                ]),
                     # output
                     (warp_fmri_wf,  datasink,  [
                                                 ("wfmri_output.warped_fmri",    "rest.{}.@warped_fmri".format(template_name)),
                                                 ("wfmri_output.wtime_filtered", "rest.{}.@time_filtered".format(template_name)),
                                                 ("wfmri_output.smooth",         "rest.{}.@smooth".format(template_name)),
                                                 ("wfmri_output.wavg_epi",       "rest.{}.@avg_epi".format(template_name)),
                                                 ("wfmri_output.warp_field",     "rest.{}.@warp_field".format(template_name)),
                                                ]),
                    ])

    if not do_group_template:
        main_wf.connect([
                         (anat_wf, warp_fmri_wf, [("anat_output.anat_biascorr",  "wfmri_input.reference_file"),
                                                  ("anat_output.warp_forward",   "wfmri_input.anat_to_mni_warp"),
                                                 ]),
                        # output
                        (warp_fmri_wf,  datasink,  [
                                                    ("wfmri_output.coreg_avg_epi", "rest.@coreg_fmri_anat"),
                                                    ("wfmri_output.coreg_others",  "rest.@coreg_others"),
                                                   ]),
                        ])

    if do_atlas:
        main_wf.connect([(anat_wf,      warp_fmri_wf, [("anat_output.atlas_anat",   "wfmri_input.atlas_anat")]),
                         (warp_fmri_wf, datasink,     [("wfmri_output.atlas_fmri",  "rest.@atlas")]),
                         ])
    return main_wf
