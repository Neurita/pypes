# -*- coding: utf-8 -*-
"""
PET-only image registration nipype workflow.
"""
import os.path as op

import nipype.pipeline.engine    as pe
from   nipype.interfaces         import io, spm
from   nipype.algorithms.misc    import Gunzip
from   nipype.interfaces.utility import IdentityInterface, Function

from   ..preproc import spm_group_template, get_bounding_box
from   ..config  import setup_node
from   .._utils  import format_pair_list
from   ..utils   import (get_datasink,
                         extend_trait_list,
                         get_input_node,
                         remove_ext,
                         get_input_file_name,
                         extension_duplicates)


def attach_spm_pet_grouptemplate(main_wf, wf_name="spm_pet_template"):
    """ Attach a PET pre-processing workflow that uses SPM12 to `main_wf`.
    This workflow picks all

    This will also attach the anat preprocessing workflow to `main_wf`. The reason
    for this is that the PET pre-processing steps here make use of anatomical MR
    pre-processing outputs.

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
    template_wf   = spm_group_template(wf_name)

    get_bbox = setup_node(Function(function=get_bounding_box,
                                   input_names=["in_file"],
                                   output_names=["bbox"]),
                          name="get_bbox")

    # the list of the raw pet subjects
    warped_pets = pe.JoinNode(interface=IdentityInterface(fields=["warped_pets"]),
                              joinsource="infosrc",
                              joinfield="warped_pets",
                              name="warped_pets")

    # warp each subject to the group template
    gunzip_template  = setup_node(Gunzip(), name="gunzip_template",)
    gunzip_pet       = setup_node(Gunzip(), name="gunzip_pet",)

    warp2template = setup_node(spm.Normalize(jobtype="estwrite", out_prefix="wgrptemplate_"),
                               name="warp2template")

    # group dataSink output substitutions
    regexp_subst = [
                     (r"/wgrptemplate{pet}_merged_mean_smooth.nii$",  "/{pet}_mni_grouptemplate.nii"),
                     (r"/w{pet}_merged_mean_smooth.nii$",             "/{pet}_mni_grouptemplate.nii"),
                   ]
    regexp_subst = format_pair_list(regexp_subst, pet=pet_fbasename)
    regexp_subst += extension_duplicates(regexp_subst)
    grp_datasink.inputs.regexp_substitutions = extend_trait_list(grp_datasink.inputs.regexp_substitutions,
                                                                 regexp_subst)

    # per-subject datasink output substitutions
    regexp_subst = [
                     (r"/{pet}_sn.mat$",           "/{pet}_warp2template_params.mat"),
                     (r"/wgrptemplate_{pet}.nii$", "/{pet}_warped2template.nii"),
                   ]
    regexp_subst = format_pair_list(regexp_subst, pet=pet_fbasename)
    regexp_subst += extension_duplicates(regexp_subst)
    datasink.inputs.regexp_substitutions = extend_trait_list(datasink.inputs.regexp_substitutions,
                                                             regexp_subst)

    # Connect the nodes
    main_wf.connect([
                     # warped pets file list input
                     (pet_wf,       warped_pets, [("pet_output.warped_files",   "warped_pets")]),

                     # group template wf
                     (warped_pets,  template_wf, [("warped_pets",               "grptemplate_input.in_files")]),

                     # template output
                     (template_wf, grp_datasink, [("grptemplate_output.template", "@pet_group_template")]),

                     # get template bounding box to apply to results
                     (template_wf, get_bbox,     [("grptemplate_output.template", "in_file")]),

                     # warp subjects to group template
                     (in_files,    gunzip_pet,      [("pet",                          "in_file")]),
                     (template_wf, gunzip_template, [("grptemplate_output.template",  "in_file")]),

                     (gunzip_pet,      warp2template, [("out_file",   "source" )]),
                     (gunzip_template, warp2template, [("out_file",   "template")]),
                     (get_bbox,        warp2template, [("bbox",       "write_bounding_box")]),

                     # output
                     (warp2template, datasink, [("normalization_parameters", "pet.@grouptemplate_warpparams"),
                                                ("normalized_source",        "pet.@grouptemplate_warped"),
                                               ]),
                   ])

    return main_wf
