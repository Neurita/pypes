# -*- coding: utf-8 -*-
"""
PET-only image registration nipype workflow.
"""
import os.path as op

import nipype.pipeline.engine    as pe
from   nipype.interfaces         import io, spm
from   nipype.algorithms.misc    import Gunzip
from   nipype.interfaces.utility import IdentityInterface, Function

from   .mrpet import attach_spm_mrpet_preprocessing
from   ..preproc import spm_group_template, get_bounding_box
from   ..config  import setup_node, get_config_setting
from   .._utils  import format_pair_list
from   ..utils   import (get_datasink,
                         extend_trait_list,
                         get_input_node,
                         remove_ext,
                         get_input_file_name,
                         extension_duplicates)


def spm_pet_grouptemplate_preproc(wf_name="spm_pet_grouptemplate_preproc"):
    """ Run the PET pre-processing workflow against the gunzip_pet.in_file files.

    This will not apply PVC or intensity normalization to the PET images.
    This will only register the PET subjects to a given custom template through the nipype input: `pet_input.pet_template`.

    For now this does not do atlas registration.

    It does:
    - SPM12 Warp PET to the given template

    Parameters
    ----------
    wf_name: str
        Name of the workflow.

    Nipype Inputs
    -------------
    pet_input.in_file: traits.File
        The raw NIFTI_GZ PET image file

    pet_input.pet_template: list of traits.File
        The template file for inter-subject registration reference.

    Nipype outputs
    --------------
    pet_output.pet_warped: existing file
        PET image normalized to the group template.

    pet_output.warp_field: existing files
        Spatial normalization parameters .mat files

    Returns
    -------
    wf: nipype Workflow
    """
    # specify input and output fields
    in_fields  = ["in_file",
                  "pet_template",]

    out_fields = ["pet_warped",
                  "warp_field",]

    # input
    pet_input = setup_node(IdentityInterface(fields=in_fields, mandatory_inputs=True),
                           name="pet_input")

    # warp each subject to the group template
    gunzip_template = setup_node(Gunzip(), name="gunzip_template",)
    gunzip_pet      = setup_node(Gunzip(), name="gunzip_pet",)

    warp2template = setup_node(spm.Normalize(jobtype="estwrite", out_prefix="wgrptemplate_"),
                               name="warp2template")

    get_bbox = setup_node(Function(function=get_bounding_box,
                                   input_names=["in_file"],
                                   output_names=["bbox"]),
                          name="get_bbox")

    # output
    pet_output = setup_node(IdentityInterface(fields=out_fields), name="pet_output")

    # Create the workflow object
    wf = pe.Workflow(name=wf_name)

    wf.connect([
                # get template bounding box to apply to results
                (pet_input,     get_bbox,        [("pet_template", "in_file")]),

                # gunzip some inputs
                (pet_input,     gunzip_pet,      [("in_file",      "in_file")]),
                (pet_input,     gunzip_template, [("pet_template", "in_file")]),

                # prepare the target parameters of the warp to template
                (gunzip_template, warp2template, [("out_file", "template")]),
                (get_bbox,        warp2template, [("bbox",     "write_bounding_box")]),

                # directly warp pet to the template
                (gunzip_pet,      warp2template, [("out_file", "source")]),

                # output
                (warp2template, pet_output, [("normalization_parameters", "warp_field"),
                                             ("normalized_files" ,        "pvc_warped"),
                                             ("normalized_source",        "pet_warped"),
                                            ]),
               ])

    return wf


def _attach_spm_pet_grouptemplate_preprocessing(main_wf, wf_name="spm_pet_grouptemplate_preproc",):
    """ Attach a PET pre-processing workflow that uses SPM12 to `main_wf`.
    This workflow picks all spm_pet_preproc outputs 'pet_output.warped_files' in `main_wf`
    to create a group template.
    This pipeline does not apply PETPVC and intensity normalization to the PET images, only PET to group
    template registration.

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

    Nipype Workflow Dependencies
    ----------------------------
    This workflow depends on:
    - spm_pet_preproc

    Returns
    -------
    main_wf: nipype Workflow
    """
    # Dependency workflows
    in_files = get_input_node(main_wf)
    datasink = get_datasink  (main_wf)

    # The base name of the 'pet' file for the substitutions
    pet_fbasename  = remove_ext(op.basename(get_input_file_name(in_files, 'pet')))

    # get the PET preprocessing pipeline
    pet_wf = spm_pet_grouptemplate_preproc(wf_name=wf_name)

    # dataSink output substitutions
    regexp_subst = [
                     (r"/w{pet}.nii", "/{pet}_grptemplate.nii"),
                   ]
    regexp_subst = format_pair_list(regexp_subst, pet=pet_fbasename)
    regexp_subst += extension_duplicates(regexp_subst)
    datasink.inputs.regexp_substitutions = extend_trait_list(datasink.inputs.regexp_substitutions,
                                                             regexp_subst)

    # Connect the nodes
    main_wf.connect([
                     # pet file input
                     (in_files, pet_wf, [("pet", "pet_input.in_file"),]),

                     (pet_wf, datasink, [("pet_output.pet_warped", "mrpet.@pet_warped"),
                                         ("pet_output.warp_field", "mrpet.@warp_field"),
                                        ]),
                     ])

    return main_wf


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
    - spm_anat_preproc if `spm_pet_template.do_petpvc` is True.

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
                                       base_directory=base_outdir,), name="group_datasink")
    grp_datasink.inputs.container = '{}_grouptemplate'.format(pet_fbasename)

    # the group template workflow
    template_wf = spm_group_template(wf_name)

    # the list of the raw pet subjects
    warped_pets = pe.JoinNode(interface=IdentityInterface(fields=["warped_pets"]),
                              joinsource="infosrc",
                              joinfield="warped_pets",
                              name="warped_pets")

    # output node
    output = setup_node(IdentityInterface(fields=["pet_template"]), name="group_template")

    # group dataSink output substitutions
    regexp_subst = [
                     (r"/wgrptemplate{pet}_merged_mean_smooth.nii$",  "/{pet}_grouptemplate_mni.nii"),
                     (r"/w{pet}_merged_mean_smooth.nii$",             "/{pet}_grouptemplate_mni.nii"),
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
                     (template_wf, output,       [("grptemplate_output.template", "pet_template")]),

                     # template output
                     (output,      grp_datasink, [("pet_template",                "@pet_group_template")]),
                   ])

    # Now we start with the correction and registration of each subject to the group template
    do_petpvc = get_config_setting('spm_pet_template.do_petpvc')

    if do_petpvc and main_wf.get_node('spm_anat_preproc') is not None:
        preproc_wf_name = "spm_mrpet_grouptemplate_preproc"
        main_wf = attach_spm_mrpet_preprocessing(main_wf, wf_name=preproc_wf_name, do_group_template=True)
    else:
        preproc_wf_name = "spm_pet_grouptemplate_preproc"
        main_wf = _attach_spm_pet_grouptemplate_preprocessing(main_wf, wf_name=preproc_wf_name)

    # add the pet template to the preproc workflow
    preproc_wf = main_wf.get_node(preproc_wf_name)
    main_wf.connect([(output, preproc_wf, [("pet_template", "pet_input.pet_template".format(preproc_wf_name))]), ])

    # per-subject datasink output substitutions
    regexp_subst = [
                     (r"/{pet}_sn.mat$",           "/{pet}_grptemplate_params.mat"),
                     (r"/wgrptemplate_{pet}.nii$", "/{pet}_grptemplate.nii"),
                   ]
    regexp_subst = format_pair_list(regexp_subst, pet=pet_fbasename)
    regexp_subst += extension_duplicates(regexp_subst)
    datasink.inputs.regexp_substitutions = extend_trait_list(datasink.inputs.regexp_substitutions,
                                                             regexp_subst)

    return main_wf

