# -*- coding: utf-8 -*-
"""
PET-MR image preprocessing nipype workflows.
"""
import os.path as op

import nipype.pipeline.engine    as pe
from   nipype.algorithms.misc    import Gunzip
from   nipype.interfaces.utility import Merge, IdentityInterface, Function
from   nipype.interfaces         import spm

from   .pvc      import petpvc_workflow
from   ..config  import (setup_node,
                        check_atlas_file,
                        get_config_setting)

from   ..preproc import (spm_normalize,
                         spm_coregister,
                         spm_apply_deformations,
                         get_bounding_box)

from   ..utils import (get_datasink,
                       spm_tpm_priors_path,
                       extend_trait_list,
                       get_input_node,
                       get_interface_node,
                       remove_ext,
                       get_input_file_name,
                       extension_duplicates)

from   .._utils import format_pair_list, flatten_list


# TODO: merge the two workflows below, maybe splitting them in
# two wf steps: pre-processing then registration.
def spm_mrpet_preprocessing(wf_name="spm_mrpet_preproc"):
    """ Run the PET pre-processing workflow against the
    gunzip_pet.in_file files.
    It depends on the anat_preproc_workflow, so if this
    has not been run, this function will run it too.

    # TODO: organize the anat2pet hack/condition somehow:
    If anat2pet:
    - SPM12 Coregister T1 and tissues to PET
    - PVC the PET image in PET space
    - SPM12 Warp PET to MNI
    else:
    - SPM12 Coregister PET to T1
    - PVC the PET image in anatomical space
    - SPM12 Warp PET in anatomical space to MNI through the
    `anat_to_mni_warp`.

    Parameters
    ----------
    wf_name: str
        Name of the workflow.

    Nipype Inputs
    -------------
    pet_input.in_file: traits.File
        The raw NIFTI_GZ PET image file

    pet_input.anat: traits.File
        Path to the high-contrast anatomical image.
        Reference file of the warp_field, i.e., the
        anatomical image in its native space.

    pet_input.anat_to_mni_warp: traits.File
        The warp field from the transformation of the
        anatomical image to the standard MNI space.

    pet_input.atlas_anat: traits.File
        The atlas file in anatomical space.

    pet_input.tissues: list of traits.File
        List of tissues files from the New Segment process.
        At least the first 3 tissues must be present.

    Nipype outputs
    --------------
    pet_output.pvc_out: existing file
        The results of the PVC process

    pet_output.brain_mask: existing file
        A brain mask calculated with the tissues file.

    pet_output.coreg_ref: existing file
        The coregistered reference image to PET space.

    pet_output.coreg_others: list of existing files
        List of coregistered files from coreg_pet.apply_to_files

    pet_output.pvc_warped: existing file
        Results from PETPVC normalized to MNI.
        The result of every internal pre-processing step
        is normalized to MNI here.

    pet_output.warp_field: existing files
        Spatial normalization parameters .mat files

    pet_output.gm_norm: existing file
        The output of the grey matter intensity
        normalization process.
        This is the last step in the PET signal correction,
        before registration.

    pet_output.atlas_pet: existing file
        Atlas image warped to PET space.
        If the `atlas_file` option is an existing file and
        `normalize_atlas` is True.

    Returns
    -------
    wf: nipype Workflow
    """
    # specify input and output fields
    in_fields  = ["in_file",
                  "anat",
                  "anat_to_mni_warp",
                  "tissues",]

    out_fields = ["brain_mask",
                  "coreg_others",
                  "coreg_ref",
                  "pvc_warped",
                  "pet_warped", # 'pet_warped' is a dummy entry to keep the fields pattern.
                  "warp_field",
                  "pvc_out",
                  "pvc_mask",
                  "gm_norm",]

    do_atlas, _ = check_atlas_file()
    if do_atlas:
        in_fields  += ["atlas_anat"]
        out_fields += ["atlas_pet" ]

    # input
    pet_input = setup_node(IdentityInterface(fields=in_fields, mandatory_inputs=True),
                           name="pet_input")

    # workflow to perform partial volume correction
    petpvc    = petpvc_workflow(wf_name="petpvc")

    merge_list = setup_node(Merge(4), name='merge_for_unzip')
    gunzipper = pe.MapNode(Gunzip(), name="gunzip", iterfield=['in_file'])

    warp_pet = setup_node(spm_normalize(), name="warp_pet")

    tpm_bbox = setup_node(Function(function=get_bounding_box,
                                   input_names=["in_file"],
                                   output_names=["bbox"]),
                          name="tpm_bbox")
    tpm_bbox.inputs.in_file = spm_tpm_priors_path()

    # output
    pet_output = setup_node(IdentityInterface(fields=out_fields), name="pet_output")

    # Create the workflow object
    wf = pe.Workflow(name=wf_name)

    # check how to perform the registration, to decide how to build the pipeline
    anat2pet = get_config_setting('registration.anat2pet', False)
    if anat2pet:
        wf.connect([
                    # inputs
                    (pet_input, petpvc,     [("in_file", "pvc_input.in_file"),
                                             ("anat",    "pvc_input.reference_file"),
                                             ("tissues", "pvc_input.tissues")]),

                    # gunzip some files for SPM Normalize
                    (petpvc,    merge_list, [("pvc_output.pvc_out",    "in1"),
                                             ("pvc_output.brain_mask", "in2"),
                                             ("pvc_output.gm_norm",    "in3")]),
                    (pet_input, merge_list, [("in_file",               "in4")]),

                    (merge_list, gunzipper, [("out", "in_file")]),

                    # warp the PET PVCed to MNI
                    (petpvc,    warp_pet,   [("pvc_output.coreg_ref", "image_to_align")]),
                    (gunzipper, warp_pet,   [("out_file",             "apply_to_files")]),
                    (tpm_bbox,  warp_pet,   [("bbox",                 "write_bounding_box")]),

                    # output
                    (petpvc,    pet_output, [("pvc_output.pvc_out",      "pvc_out"),
                                             ("pvc_output.brain_mask",   "brain_mask"),
                                             ("pvc_output.coreg_ref",    "coreg_ref"),
                                             ("pvc_output.coreg_others", "coreg_others"),
                                             ("pvc_output.gm_norm",      "gm_norm")]),

                    # output
                    (warp_pet,  pet_output, [("normalized_files",  "pvc_warped"),
                                             ("deformation_field", "warp_field")]),
                   ])
    else: # PET 2 ANAT
        collector  = setup_node(Merge(2), name='merge_for_warp')
        apply_warp = setup_node(spm_apply_deformations(), name="warp_pet")

        wf.connect([
                    # inputs
                    (pet_input, petpvc,     [("in_file", "pvc_input.in_file"),
                                             ("anat",    "pvc_input.reference_file"),
                                             ("tissues", "pvc_input.tissues")]),

                    # gunzip some files for SPM Normalize
                    (petpvc,    merge_list, [("pvc_output.pvc_out",    "in1"),
                                             ("pvc_output.brain_mask", "in2"),
                                             ("pvc_output.gm_norm",    "in3")]),
                    (pet_input, merge_list, [("in_file",               "in4")]),

                    (merge_list, gunzipper, [("out",                   "in_file")]),

                    # warp the PET PVCed to MNI
                    (gunzipper,   collector,   [("out_file",             "in1")]),
                    (petpvc,      collector,   [("pvc_output.coreg_ref", "in2")]),

                    (pet_input,   apply_warp,  [("anat_to_mni_warp", "deformation_file")]),
                    (collector,   apply_warp,  [("out",              "apply_to_files")]),
                    (tpm_bbox,    apply_warp,  [("bbox",             "write_bounding_box")]),

                    # output
                    (petpvc,    pet_output, [("pvc_output.pvc_out",      "pvc_out"),
                                             ("pvc_output.brain_mask",   "brain_mask"),
                                             ("pvc_output.petpvc_mask",  "petpvc_mask"),
                                             ("pvc_output.coreg_ref",    "coreg_ref"),
                                             ("pvc_output.coreg_others", "coreg_others"),
                                             ("pvc_output.gm_norm",      "gm_norm")]),

                    # output
                    (apply_warp,  pet_output, [("normalized_files",  "pvc_warped"),
                                               ("deformation_field", "warp_field")]),
                   ])


    if do_atlas:
        coreg_atlas = setup_node(spm_coregister(cost_function="mi"), name="coreg_atlas")

        # set the registration interpolation to nearest neighbour.
        coreg_atlas.inputs.write_interp = 0
        wf.connect([
                    (pet_input,   coreg_atlas, [("anat",                 "source")]),
                    (petpvc,      coreg_atlas, [("pvc_output.coreg_ref", "target")]),
                    (pet_input,   coreg_atlas, [("atlas_anat",           "apply_to_files")]),
                    (coreg_atlas, pet_output,  [("coregistered_files",   "atlas_pet")]),
        ])

    return wf


def spm_mrpet_grouptemplate_preprocessing(wf_name="spm_mrpet_grouptemplate_preproc"):
    """ Run the PET pre-processing workflow against the gunzip_pet.in_file files.
    It depends on the anat_preproc_workflow, so if this has not been run, this function
    will run it too.

    This is identical to the workflow defined in `spm_mrpet_preprocessing`,
    with the only difference that we now normalize all subjects agains a custom
    template using the spm Old Normalize interface.

    It does:
    - SPM12 Coregister T1 and tissues to PET
    - PVC the PET image in PET space
    - SPM12 Warp PET to the given template

    Parameters
    ----------
    wf_name: str
        Name of the workflow.

    Nipype Inputs
    -------------
    pet_input.in_file: traits.File
        The raw NIFTI_GZ PET image file.

    pet_input.atlas_anat: traits.File
        The atlas file in anatomical space.

    pet_input.anat: traits.File
        Path to the high-contrast anatomical image.
        Reference file of the warp_field, i.e., the anatomical image in its native space.

    pet_input.tissues: list of traits.File
        List of tissues files from the New Segment process. At least the first
        3 tissues must be present.

    pet_input.pet_template: traits.File
        The template file for inter-subject registration reference.

    Nipype outputs
    --------------
    pet_output.pvc_out: existing file
        The results of the PVC process.

    pet_output.brain_mask: existing file
        A brain mask calculated with the tissues file.

    pet_output.coreg_ref: existing file
        The coregistered reference image to PET space.

    pet_output.coreg_others: list of existing files
        List of coregistered files from coreg_pet.apply_to_files.

    pet_output.pet_warped: existing file
        PET image normalized to the group template.

    pet_output.pvc_warped: existing file
        The outputs of the PETPVC workflow normalized to the group template.
        The result of every internal pre-processing step is normalized to the
        group template here.

    pet_output.warp_field: existing files
        Spatial normalization parameters .mat files.

    pet_output.gm_norm: existing file
        The output of the grey matter intensity normalization process.
        This is the last step in the PET signal correction, before registration.

    pet_output.atlas_pet: existing file
        Atlas image warped to PET space.
        If the `atlas_file` option is an existing file and `normalize_atlas` is True.

    Returns
    -------
    wf: nipype Workflow
    """
    # specify input and output fields
    in_fields  = ["in_file",
                  "anat",
                  "tissues",
                  "pet_template"]

    out_fields = ["brain_mask",
                  "coreg_others",
                  "coreg_ref",
                  "pvc_warped",
                  "pet_warped",
                  "warp_field",
                  "pvc_out",
                  "pvc_mask",
                  "gm_norm",]

    do_atlas, _ = check_atlas_file()
    if do_atlas:
        in_fields  += ["atlas_anat"]
        out_fields += ["atlas_pet" ]

    # input
    pet_input = setup_node(IdentityInterface(fields=in_fields, mandatory_inputs=True),
                           name="pet_input")

    # workflow to perform partial volume correction
    petpvc = petpvc_workflow(wf_name="petpvc")

    unzip_mrg = setup_node(Merge(4), name='merge_for_unzip')
    gunzipper = pe.MapNode(Gunzip(), name="gunzip", iterfield=['in_file'])

    # warp each subject to the group template
    gunzip_template = setup_node(Gunzip(), name="gunzip_template",)
    gunzip_pet      = setup_node(Gunzip(), name="gunzip_pet",)

    warp_mrg = setup_node(Merge(2), name='merge_for_warp')
    warp2template = setup_node(spm.Normalize(jobtype="estwrite", out_prefix="wgrptemplate_"),
                               name="warp2template",)

    get_bbox = setup_node(Function(function=get_bounding_box,
                                   input_names=["in_file"],
                                   output_names=["bbox"]),
                          name="get_bbox")

    # output
    pet_output = setup_node(IdentityInterface(fields=out_fields), name="pet_output")

    # Create the workflow object
    wf = pe.Workflow(name=wf_name)

    wf.connect([
                # inputs
                (pet_input,   petpvc,  [("in_file", "pvc_input.in_file"),
                                        ("anat",    "pvc_input.reference_file"),
                                        ("tissues", "pvc_input.tissues")]),

                # get template bounding box to apply to results
                (pet_input, get_bbox,  [("pet_template", "in_file")]),

                # gunzip some inputs
                (pet_input, gunzip_pet,      [("in_file",      "in_file")]),
                (pet_input, gunzip_template, [("pet_template", "in_file")]),

                # gunzip some files for SPM Normalize
                (petpvc,    unzip_mrg, [("pvc_output.pvc_out",    "in1"),
                                        ("pvc_output.brain_mask", "in2"),
                                        ("pvc_output.gm_norm",    "in3")]),
                (pet_input, unzip_mrg, [("in_file",               "in4")]),

                (unzip_mrg, gunzipper, [("out", "in_file")]),

                (gunzipper, warp_mrg,  [("out_file", "in1")]),

                (warp_mrg, warp2template, [(("out", flatten_list), "apply_to_files")]),

                # prepare the target parameters of the warp to template
                (gunzip_pet,      warp2template, [("out_file", "source")]),
                (gunzip_template, warp2template, [("out_file", "template")]),
                (get_bbox,        warp2template, [("bbox",     "write_bounding_box")]),

                # output
                (warp2template, pet_output, [("normalization_parameters", "warp_field"),
                                             ("normalized_files" ,        "pvc_warped"),
                                             ("normalized_source",        "pet_warped"),
                                            ]),

                # output
                (petpvc,   pet_output, [("pvc_output.pvc_out",      "pvc_out"),
                                        ("pvc_output.brain_mask",   "brain_mask"),
                                        ("pvc_output.coreg_ref",    "coreg_ref"),
                                        ("pvc_output.coreg_others", "coreg_others"),
                                        ("pvc_output.gm_norm",      "gm_norm")]),
                ])

    if do_atlas:
        coreg_atlas = setup_node(spm_coregister(cost_function="mi"), name="coreg_atlas")

        # set the registration interpolation to nearest neighbour.
        coreg_atlas.inputs.write_interp = 0
        wf.connect([
                    (pet_input,   coreg_atlas, [("anat",                 "source")]),
                    (petpvc,      coreg_atlas, [("pvc_output.coreg_ref", "target")]),
                    (pet_input,   coreg_atlas, [("atlas_anat",           "apply_to_files")]),
                    (coreg_atlas, pet_output,  [("coregistered_files",   "atlas_pet")]),

                    # warp the atlas to the template space as well
                    (coreg_atlas, warp_mrg,    [("coregistered_files",   "in2")]),
        ])

    return wf


def attach_spm_mrpet_preprocessing(main_wf, wf_name="spm_mrpet_preproc",
                                   do_group_template=False):
    """ Attach a PET pre-processing workflow that uses SPM12 to `main_wf`.
    This workflow needs MRI based

    Nipype Inputs for `main_wf`
    ---------------------------
    Note: The `main_wf` workflow is expected to have an `input_files` and a
    `datasink` nodes.

    input_files.select.pet: input node

    datasink: nipype Node

    Parameters
    ----------
    main_wf: nipype Workflow

    wf_name: str
        Name of the preprocessing workflow

    do_group_template: bool
        If True will attach the group template creation and pre-processing pipeline.

    Nipype Workflow Dependencies
    ----------------------------
    This workflow depends on:
    - spm_anat_preproc

    Returns
    -------
    main_wf: nipype Workflow
    """
    # Dependency workflows
    in_files = get_input_node(main_wf)
    datasink = get_datasink  (main_wf)

    anat_output = get_interface_node(main_wf, "anat_output")

    # The base name of the 'pet' file for the substitutions
    anat_fbasename = remove_ext(op.basename(get_input_file_name(in_files, 'anat')))
    pet_fbasename  = remove_ext(op.basename(get_input_file_name(in_files, 'pet')))

    # get the PET preprocessing pipeline
    if do_group_template:
        pet_wf = spm_mrpet_grouptemplate_preprocessing(wf_name=wf_name)
        template_name = 'grptemplate'
        output_subfolder = 'grp_template'
    else:
        pet_wf = spm_mrpet_preprocessing(wf_name=wf_name)
        template_name = 'mni'
        output_subfolder = 'std_template'

    # dataSink output substitutions
    regexp_subst = [
                     (r"/{pet}_.*_pvc.nii.gz$",           "/{pet}_pvc.nii.gz"),
                     (r"/{pet}_.*_pvc_maths.nii.gz$",     "/{pet}_pvc_norm.nii.gz"),
                     (r"/{pet}_.*_pvc_intnormed.nii.gz$", "/{pet}_pvc_norm.nii.gz"),
                     (r"/tissues_brain_mask.nii$",        "/brain_mask_anat.nii"),
                     (r"/w{pet}.nii",                     "/{pet}_{template}.nii"),
                     (r"/w{pet}_.*_pvc.nii$",             "/{pet}_pvc_{template}.nii"),
                     (r"/w{pet}_.*_pvc_maths.nii$",       "/{pet}_pvc_norm_{template}.nii"),
                     (r"/w{pet}_.*_pvc_intnormed.nii$",   "/{pet}_pvc_norm_{template}.nii"),
                     (r"/wbrain_mask.nii",                "/brain_mask_{template}.nii"),
                     (r"/r{pet}.nii",                     "/{pet}_anat.nii"),
                     (r"/r{pet}_.*_pvc.nii$",             "/{pet}_pvc_anat.nii"),
                     (r"/r{pet}_.*_pvc_maths.nii$",       "/{pet}_pvc_norm_anat.nii"),
                     (r"/r{pet}_.*_pvc_intnormed.nii$",   "/{pet}_pvc_norm_anat.nii"),
                     (r"/y_rm{anat}_corrected.nii",       "/{anat}_{pet}_warpfield.nii"),
                     (r"/rm{anat}_corrected.nii$",        "/{anat}_{pet}.nii"),
                     (r"/rc1{anat}_corrected.nii$",       "/gm_{pet}.nii"),
                     (r"/rc2{anat}_corrected.nii$",       "/wm_{pet}.nii"),
                     (r"/rc3{anat}_corrected.nii$",       "/csf_{pet}.nii"),
                   ]
    regexp_subst = format_pair_list(regexp_subst, pet=pet_fbasename, anat=anat_fbasename,
                                    template=template_name)

    # prepare substitution for atlas_file, if any
    do_atlas, atlas_file = check_atlas_file()
    if do_atlas:
        atlas_basename = remove_ext(op.basename(atlas_file))
        regexp_subst.extend([
                             (r"/[\w]*{atlas}\.nii$", "/{atlas}_{pet}.nii"),
                            ])
        regexp_subst = format_pair_list(regexp_subst, pet=pet_fbasename, atlas=atlas_basename)

    regexp_subst += extension_duplicates(regexp_subst)

    datasink.inputs.regexp_substitutions = extend_trait_list(datasink.inputs.regexp_substitutions,
                                                             regexp_subst)

    # Connect the nodes
    main_wf.connect([
                     # pet file input
                     (in_files, pet_wf, [("pet", "pet_input.in_file")]),

                     # pet to anat registration
                     (anat_output,  pet_wf, [("anat_biascorr",  "pet_input.anat"),
                                             ("tissues_native", "pet_input.tissues"),
                                            ]),

                     (pet_wf, datasink, [
                                         ("pet_output.gm_norm",      "mrpet.@norm"),
                                         ("pet_output.coreg_others", "mrpet.tissues"),
                                         ("pet_output.coreg_ref",    "mrpet.@anat"),
                                         ("pet_output.pvc_mask",     "mrpet.@pvc_mask"),
                                         ("pet_output.pvc_out",      "mrpet.@pvc"),
                                         ("pet_output.brain_mask",   "mrpet.@brain_mask"),
                                         ("pet_output.pvc_warped",   "mrpet.{}.@pvc".format(output_subfolder)),
                                         ("pet_output.warp_field",   "mrpet.{}.@warp_field".format(output_subfolder)),
                                         ("pet_output.pet_warped",   "mrpet.{}.@pet_warped".format(output_subfolder)),
                                        ]),
                     ])

    if not do_group_template:
        # Connect the nodes
        main_wf.connect([
                         # pet to anat registration
                         (anat_output,  pet_wf, [("warp_forward", "pet_input.anat_to_mni_warp"),]),
                        ])


    if do_atlas:
            main_wf.connect([(anat_output, pet_wf,   [("atlas_anat",           "pet_input.atlas_anat")]),
                             (pet_wf,      datasink, [("pet_output.atlas_pet", "mrpet.@atlas")]),
                            ])

    return main_wf
