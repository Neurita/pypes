# -*- coding: utf-8 -*-
"""
PET group template registration nipype workflow.
"""
import os.path as op

import nipype.pipeline.engine as pe
from   nipype.interfaces.io import DataSink
from   nipype.interfaces import IdentityInterface

from   .mrpet import attach_spm_mrpet_preprocessing
from   ..preproc import (spm_create_group_template_wf,
                         spm_register_to_template_wf,)
from   ..config  import setup_node, get_config_setting
from   .._utils  import format_pair_list, flatten_list
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
    grp_datasink = pe.Node(DataSink(parameterization=False,
                                    base_directory=base_outdir,),
                                    name='{}_grouptemplate_datasink'.format(pet_fbasename))
    grp_datasink.inputs.container = '{}_grouptemplate'.format(pet_fbasename)

    # the list of the raw pet subjects
    warped_pets = pe.JoinNode(interface=IdentityInterface(fields=["warped_pets"]),
                              joinsource="infosrc",
                              joinfield="warped_pets",
                              name="warped_pets")

    # the group template workflow
    template_wf = spm_create_group_template_wf(wf_name)

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
                     (pet_wf,       warped_pets, [("warp_output.warped_files",    "warped_pets")]),

                     # group template wf
                     (warped_pets,  template_wf, [(("warped_pets", flatten_list), "grptemplate_input.in_files")]),

                     # output node
                     (template_wf, output,       [("grptemplate_output.template", "pet_template")]),

                     # template output
                     (output,      grp_datasink, [("pet_template",                "@pet_group_template")]),
                   ])

    # Now we start with the correction and registration of each subject to the group template
    do_petpvc = get_config_setting('spm_pet_template.do_petpvc')
    if do_petpvc:
        if main_wf.get_node('spm_anat_preproc') is None:
            raise AttributeError("Expected `spm_anat_preproc` workflow node to attach PETPVC.")

        preproc_wf_name = "spm_mrpet_grouptemplate_preproc"
        main_wf = attach_spm_mrpet_preprocessing(main_wf, wf_name=preproc_wf_name, do_group_template=True)

        preproc_wf = main_wf.get_node(preproc_wf_name)
        main_wf.connect([(output, preproc_wf, [("pet_template", "pet_input.pet_template".format(preproc_wf_name))]), ])
    else:
        # add the pet template to the preproc workflow
        reg_wf = spm_register_to_template_wf(wf_name="spm_pet_register_to_grouptemplate")
        main_wf.connect([
                         (output,   reg_wf,  [("pet_template",  "reg_input.template")]),
                         (in_files, reg_wf,  [("pet",           "reg_input.in_file"),]),

                         (reg_wf,   datasink, [("reg_output.warped",     "pet.grp_template.@warped"),
                                               ("reg_output.warp_field", "pet.grp_template.@warp_field"),
                                              ]),
                         ])

    # per-subject datasink output substitutions
    regexp_subst = [
                     (r"/{pet}_sn.mat$",           "/{pet}_grptemplate_params.mat"),
                     (r"/wgrptemplate_{pet}.nii$", "/{pet}_grptemplate.nii"),
                     (r"/w{pet}.nii",              "/{pet}_grptemplate.nii"),
                   ]
    regexp_subst = format_pair_list(regexp_subst, pet=pet_fbasename)
    regexp_subst += extension_duplicates(regexp_subst)
    datasink.inputs.regexp_substitutions = extend_trait_list(datasink.inputs.regexp_substitutions,
                                                             regexp_subst)

    return main_wf

