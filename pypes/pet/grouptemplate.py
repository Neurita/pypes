# -*- coding: utf-8 -*-
"""
PET-only image registration nipype workflow.
"""
import os.path as op

import nipype.pipeline.engine    as pe
from   nipype.interfaces         import io, spm
from   nipype.algorithms.misc    import Gunzip
from   nipype.interfaces.utility import IdentityInterface, Function, Select, Merge

from   ..preproc import spm_group_template, get_bounding_box
from   ..config  import setup_node, get_config_setting
from   .._utils  import format_pair_list
from   ..utils   import (get_datasink,
                         extend_trait_list,
                         get_input_node,
                         remove_ext,
                         get_input_file_name,
                         extension_duplicates)


def attach_spm_pet_grouptemplate(main_wf, wf_name="spm_pet_template"):
    """ Attach a PET pre-processing workflow that uses SPM12 to `main_wf`.
    This workflow picks all spm_pet_preproc outputs 'pet_output.warped_files' in `main_wf`
    to create a group template.

    Parameters
    ----------
    main_wf: nipype Workflow

    wf_name: str
        Name of the preprocessing workflow

    Nipype Inputs for `main_wf`
    ---------------------------
    Note: The `main_wf` workflow is expected to have an `input_files` and a `datasink` nodes.

    pet_output.warped_files: input node

    datasink: nipype Node

    spm_pet_preproc: nipype Workflow

    Nipype Outputs
    --------------
    group_template.pet_template: file
        The path to the PET group template.

    Nipype Workflow Dependencies
    ----------------------------
    This workflow depends on:
    - spm_pet_preproc

    Returns
    -------
    main_wf: nipype Workflow
    """
    # Dependency workflows
    pet_wf = main_wf.get_node("spm_pet_preproc")

    in_files = get_input_node(main_wf)
    datasink = get_datasink(main_wf, name='datasink')

    # The base name of the 'pet' file for the substitutions
    pet_fbasename = remove_ext(op.basename(get_input_file_name(in_files, 'pet')))

    # the group template datasink
    base_outdir  = datasink.inputs.base_directory
    grp_datasink = pe.Node(io.DataSink(parameterization=False,
                                       base_directory=base_outdir,),
                           name="group_datasink")
    grp_datasink.inputs.container = '{}_grouptemplate'.format(pet_fbasename)

    # the group template workflow
    template_wf = spm_group_template(wf_name)

    # the list of the raw pet subjects
    warped_pets = pe.JoinNode(interface=IdentityInterface(fields=["warped_pets"]),
                              joinsource="infosrc",
                              joinfield="warped_pets",
                              name="warped_pets")

    # output node
    output = setup_node(IdentityInterface(fields=["pet_template"]),
                        name="group_template")

    # group dataSink output substitutions
    regexp_subst = [
                     (r"/wgrptemplate{pet}_merged_mean_smooth.nii$",  "/{pet}_mni_grouptemplate.nii"),
                     (r"/w{pet}_merged_mean_smooth.nii$",             "/{pet}_mni_grouptemplate.nii"),
                   ]
    regexp_subst = format_pair_list(regexp_subst, pet=pet_fbasename)
    regexp_subst += extension_duplicates(regexp_subst)
    grp_datasink.inputs.regexp_substitutions = extend_trait_list(grp_datasink.inputs.regexp_substitutions,
                                                                 regexp_subst)

    # Connect the nodes
    main_wf.connect([
                     # warped pets file list input
                     (pet_wf,       warped_pets, [("pet_output.warped_files",     "warped_pets")]),

                     # group template wf
                     (warped_pets,  template_wf, [("warped_pets",                 "grptemplate_input.in_files")]),

                     # output node
                     (template_wf, output,       [("grptemplate_output.template", "group_template.pet_template")]),

                     # template output
                     (output,      grp_datasink, [("group_template.pet_template", "@pet_group_template")]),
                   ])

    return main_wf


def attach_spm_pet_grouptemplate_preproc(main_wf, wf_name='spm_pet_grouptemplate_preproc'):
    """ Attach a PET pre-processing workflow that uses SPM12 to `main_wf`.
    This workflow picks all spm_pet_preproc outputs 'pet_output.warped_files' in `main_wf`
    to create a group template.

    Parameters
    ----------
    main_wf: nipype Workflow

    wf_name: str
        Name of the preprocessing workflow

    Nipype Inputs
    -------------
    datasink: nipype Node

    spm_pet_preproc: nipype Workflow

    spm_pet_template: nipype Workflow

    spm_anat_preproc: nipype Node
        If spm_pet_template.do_petpvc setting is True.
    """
    # Dependency workflows
    do_petpvc = get_config_setting('spm_pet_template.do_petpvc')
    if do_petpvc:
        anat_wf  = main_wf.get_node("spm_anat_preproc")

    in_files = get_input_node(main_wf)
    datasink = get_datasink(main_wf, name='datasink')

    # The base name of the 'pet' file for the substitutions
    pet_fbasename = remove_ext(op.basename(get_input_file_name(in_files, 'pet')))

    # warp each subject to the group template
    gunzip_template = setup_node(Gunzip(), name="gunzip_template",)
    gunzip_pet      = setup_node(Gunzip(), name="gunzip_pet",)

    warp2template = setup_node(spm.Normalize(jobtype="estwrite", out_prefix="wgrptemplate_"),
                               name="warp2template")

    get_bbox = setup_node(Function(function=get_bounding_box,
                                   input_names=["in_file"],
                                   output_names=["bbox"]),
                          name="get_bbox")

    # per-subject datasink output substitutions
    regexp_subst = [
                     (r"/{pet}_sn.mat$",           "/{pet}_warp2template_params.mat"),
                     (r"/wgrptemplate_{pet}.nii$", "/{pet}_warped2template.nii"),
                   ]
    regexp_subst = format_pair_list(regexp_subst, pet=pet_fbasename)
    regexp_subst += extension_duplicates(regexp_subst)
    datasink.inputs.regexp_substitutions = extend_trait_list(datasink.inputs.regexp_substitutions,
                                                             regexp_subst)

    main_wf.connect([
                     # get template bounding box to apply to results
                     (main_wf,     get_bbox,     [("group_template.pet_template",     "in_file")]),

                     # gunzip some inputs
                     (in_files,    gunzip_pet,      [("pet",                          "in_file")]),
                     (main_wf,     gunzip_template, [("group_template.pet_template",  "in_file")]),

                    # prepare the target parameters of the warp to template
                     (gunzip_template, warp2template, [("out_file", "template")]),
                     (get_bbox,        warp2template, [("bbox",     "write_bounding_box")]),

                     # output
                     (warp2template, datasink, [("normalization_parameters", "pet.@grouptemplate_warpparams"),
                                                ("normalized_source",        "pet.@grouptemplate_warped"),
                                               ]),
                   ])

    if not do_petpvc:
        # Connect the nodes
        main_wf.connect([
                         # directly warp pet to the template
                         (gunzip_pet,      warp2template, [("out_file", "source" )]),
                       ])

    else: # do petpvc before registration

        from   ..preproc import (spm_apply_deformations,
                                 spm_normalize,
                                 spm_coregister,
                                 petpvc_cmd,
                                 petpvc_mask,
                                 intensity_norm)

        #
        # import nipype.pipeline.engine    as pe
        # from   nipype.algorithms.misc    import Gunzip
        # from   nipype.interfaces.utility import Select, Merge, IdentityInterface
        #
        # from   ..config  import setup_node, check_atlas_file
        # from   ..preproc import (spm_apply_deformations,
        #                          spm_normalize,
        #                          spm_coregister,
        #                          petpvc_cmd,
        #                          petpvc_mask,
        #                          intensity_norm)
        #
        # from   ..utils import (get_datasink,
        #                        extend_trait_list,
        #                        get_input_node,
        #                        remove_ext,
        #                        get_input_file_name,
        #                        extension_duplicates)
        #
        # from   .._utils import (flatten_list,
        #                         format_pair_list)

        # fixed parameters of the NUK mMR
        psf_fwhm = (4.3, 4.3, 4.3)

        # coreg pet
        coreg_pet   = setup_node(spm_coregister(cost_function="mi"), name="coreg_pet")
        #tissues_sel = setup_node(Select(index=[0, 1, 2]),            name="tissues")
        #select_gm   = setup_node(Select(index=[0]),                  name="select_gm")
        #rbvpvc      = setup_node(petpvc_cmd(fwhm_mm=psf_fwhm,
        #                                    pvc_method='RBV'),       name="rbvpvc")
        #warp_pet    = setup_node(spm_apply_deformations(),           name="warp_pet")
        #merge_lists = setup_node(Merge(2),                           name='merge_for_warp')

        #unzip_mrg = setup_node(Merge(3),                             name='merge_for_unzip')
        #gunzipper = pe.MapNode(Gunzip(),                             name="gunzip", iterfield=['in_file'])

        # workflow to create the mask
        #mask_wf = petpvc_mask(wf_name="petpvc_mask")

        # workflow for intensity normalization
        #norm_wf = intensity_norm(wf_name="intensity_norm_gm")

        # Connect the nodes
        main_wf.connect([
                         (gunzip_pet,      warp2template, [("out_file",   "source" )]),
        ])

        # Connect the nodes
        main_wf.connect([
                         (pet_wf, datasink, [
                                        ("pet_output.out_file",     "mrpet.@pvc"),
                                        ("pet_output.coreg_others", "mrpet.others"),
                                        ("pet_output.coreg_pet",    "mrpet.@anat"),
                                        ("pet_output.brain_mask",   "mrpet.@brain_mask"),
                                        ("pet_output.gm_norm",      "mrpet.@norm"),
                                        ("pet_output.mni_pet",      "mrpet.warped2mni"),
                                       ]),

                         # warp the PET PVCed to MNI
                         (gunzipper,  merge_lists, [("out_file",            "in1")]),
                         (coreg_pet,  merge_lists, [("coregistered_source", "in2")]),

                         (merge_lists, warp_pet,   [("out",                 "apply_to_files")]),

                        # pet to anat registration
                        (anat_wf,  pet_wf, [("new_segment.bias_corrected_images",     "pet_input.coreg_target"),
                                            ("new_segment.native_class_images",       "pet_input.tissues"),
                                           ]),
                        ])



    wf.connect([
                # inputs
                (pet_input,   gunzip_pet,  [("in_file",            "in_file")]),
                (pet_input,   coreg_pet,   [("reference_file",     "target")]),
                (pet_input,   warp_pet,    [("warp_field",         "deformation_file")]),
                (pet_input,   tissues_sel, [("tissues",            "inlist")]),

                # unzip to coregister the PET image to anatomical image. 'coreg_pet.target' is an input for this wf.
                (gunzip_pet,  coreg_pet,  [("out_file",            "source")]),

                # the list of tissues to the mask wf and the GM for PET intensity normalization
                (tissues_sel, select_gm,  [(("out", flatten_list), "inlist")]),
                (tissues_sel, mask_wf,    [(("out", flatten_list), "split_tissues.inlist")]),

                # the coregistered PET to PVC correction
                (coreg_pet,   rbvpvc,     [("coregistered_source", "in_file")]),

                # the merged file with 4 tissues to PCV correction
                (mask_wf,     rbvpvc,     [("merge_tissues.merged_file", "mask_file")]),

                # normalize voxel values of PET PVCed by demeaning whole by GM PET voxel values
                (rbvpvc,      norm_wf,    [("out_file",            "mean_value.in_file")]),
                (rbvpvc,      norm_wf,    [("out_file",            "gm_norm.in_file")]),
                (select_gm,   norm_wf,    [("out",                 "mean_value.mask_file")]),

                # gunzip some files for SPM Normalize12
                (rbvpvc,      unzip_mrg,  [("out_file",            "in1")]),
                (mask_wf,     unzip_mrg,  [("brain_mask.out_file", "in2")]),
                (norm_wf,     unzip_mrg,  [("gm_norm.out_file",    "in3")]),
                (unzip_mrg,   gunzipper,  [("out",                 "in_file")]),

                # warp the PET PVCed to MNI
                (gunzipper,  merge_lists, [("out_file",            "in1")]),
                (coreg_pet,  merge_lists, [("coregistered_source", "in2")]),

                (merge_lists, warp_pet,   [("out",                 "apply_to_files")]),

                # output
                (rbvpvc,      pet_output, [("out_file",            "out_file")]),
                (mask_wf,     pet_output, [("brain_mask.out_file", "brain_mask")]),
                (coreg_pet,   pet_output, [("coregistered_source", "coreg_pet")]),
                (coreg_pet,   pet_output, [("coregistered_files",  "coreg_others")]),
                (norm_wf,     pet_output, [("gm_norm.out_file",    "gm_norm")]),
                (warp_pet,    pet_output, [("normalized_files",    "mni_pet")]),
                (warp_pet,    pet_output, [("deformation_field",   "warp_field")]),
               ])


    return main_wf


