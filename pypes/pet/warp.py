# -*- coding: utf-8 -*-
"""
PET-only image registration nipype workflow.
"""
import os.path as op

import nipype.pipeline.engine    as pe
from   nipype.interfaces         import fsl, io, spm
from   nipype.algorithms.misc    import Gunzip
from   nipype.interfaces.utility import Function, Select, Merge, IdentityInterface

from   ..preproc import (spm_apply_deformations,
                         spm_normalize,
                         spm_coregister,
                         petpvc_cmd,
                         petpvc_mask,
                         intensity_norm)

from   ..utils import (setup_node,
                       get_datasink,
                       extend_trait_list,
                       get_input_node,
                       remove_ext,
                       get_input_file_name,
                       extension_duplicates)

from   .._utils import (flatten_list,
                        format_pair_list)


def printit(x):
    print('====================================================')
    print(x)
    print('====================================================')
    return x


def spm_group_template(wf_name="spm_group_template"):
    """ Pick all subject files in `grptemplate_input.in_files`, calculate an average
    image and smooth it with `"{}_smooth".format(wf_name)` node (you can configure the smooth `fwhm` from
    a config file.).

    It does:
    - calculate a mean image (across subjects) and
    - smooth it with 8x8x8 gaussian kernel -> this is the template.
    - Finally, warp all PET images again to this template.
    - If tissue data is available from MR, do PVC.

    Parameters
    ----------
    wf_name: str
        Name of the workflow.

    Nipype Inputs
    -------------
    grptemplate_input.in_files: list of traits.File
        The raw NIFTI_GZ PET image files

    Nipype outputs
    --------------
    grptemplate_output.template: existing file
        The common custom PET template file.

    Returns
    -------
    wf: nipype Workflow
    """
    # input
    input = setup_node(IdentityInterface(fields=["in_files"]),
                       name="grptemplate_input",
                       )

    # merge
    merge  = setup_node(fsl.Merge(dimension='t'),
                        name='merge_time')

    # average
    average = setup_node(fsl.MeanImage(dimension='T'),
                         name='average')

    # smooth
    smooth = setup_node(fsl.IsotropicSmooth(fwhm=8),
                        name="{}_smooth".format(wf_name))

    #output
    output = setup_node(IdentityInterface(fields=["template"]),
                       name="grptemplate_output",
                       )

    # Create the workflow object
    wf = pe.Workflow(name=wf_name)

    wf.connect([
                # input
                (input,       merge,   [("in_files",         "in_files")]),

                # merge, average and smooth
                (merge,       average, [("merged_file",      "in_file")]),
                (average,     smooth,  [("out_file",         "in_file")]),

                # output
                (smooth,     output,   [("out_file", "template")]),
               ])

    return wf


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

    input_files.select.pet: input node

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

    base_outdir  = datasink.inputs.base_directory
    grp_datasink = pe.Node(io.DataSink(parameterization=False,
                                       base_directory=base_outdir,),
                           name="group_datasink")

    # The base name of the 'pet' file for the substitutions
    pet_fbasename = remove_ext(op.basename(get_input_file_name(in_files, 'pet')))

    template_wf   = spm_group_template(wf_name)

    # the list of the raw pet subjects
    warped_pets = pe.JoinNode(interface=IdentityInterface(fields=["warped_pets"]),
                              joinsource="infosrc",
                              joinfield="warped_pets",
                              name="warped_pets")

    # warp each subject to the group template
    gunzip_template  = setup_node(Gunzip(), name="gunzip_template",)
    gunzip_pet       = setup_node(Gunzip(), name="gunzip_pet",)

    warp2template = setup_node(spm.Normalize(jobtype="estwrite", out_prefix="wgrptemplate_"),
                               #type="map",
                               name="warp2template")
                               #iterfield=["source"])

    print_out = pe.Node(Function(input_names=['x'], output_names=['x'], function=printit),
                        name='print')

    # group dataSink output substitutions
    regexp_subst = [
                     (r"/wgrptemplate{pet}_merged_mean_smooth.nii$",  "/{pet}_mni_grouptemplate.nii"),
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
                     (template_wf, grp_datasink, [("grptemplate_output.template", "pet_group_template")]),

                     # warp subjects to group template
                     #(template_wf,    print_out, [("grptemplate_output.template",   "x")]),
                     (in_files,    gunzip_pet,      [("pet",                          "in_file")]),
                     (template_wf, gunzip_template, [("grptemplate_output.template",  "in_file")]),

                     (gunzip_pet,      warp2template, [("out_file",   "source" )]),
                     (gunzip_template, warp2template, [("out_file",   "template")]),

                     # output
                     (warp2template, datasink, [("normalization_parameters", "pet.@grouptemplate_warpparams"),
                                                ("normalized_source",        "pet.@grouptemplate_warped"),
                                               ]),
                   ])

    return main_wf


def spm_pet_preproc(wf_name="spm_pet_preproc"):
    """ Run a PET-only pre-processing workflow against the gunzip_pet.in_file files.
    It depends on the anat_preproc_workflow, so if this has not been run, this function
    will run it too.

    It does:
    - Warp each individual PET image to the default (SPM) PET template (H2O),

    Parameters
    ----------
    wf_name: str
        Name of the workflow.

    Nipype Inputs
    -------------
    pet_input.in_files: list of traits.File
        The raw NIFTI_GZ PET image files

    pet_input.tissues: list of traits.File



    Nipype outputs
    --------------
    pet_output.warped_files: list of existing file
        The warped PET files.

    Returns
    -------
    wf: nipype Workflow
    """
    # input
    pet_input = setup_node(IdentityInterface(fields=["in_files"]),
                                             name="pet_input",
                                             )

    gunzip_pet  = setup_node(Gunzip(), name="gunzip_pet",
                             )

    warp = setup_node(spm_normalize(voxel_size=[2, 2, 2]),
                      name="pet_normalize12")

    # output
    pet_output = setup_node(IdentityInterface(fields=["warped_files",
                                                     ]),
                                               name="pet_output")

    # Create the workflow object
    wf = pe.Workflow(name=wf_name)

    wf.connect([
                # inputs
                (pet_input,   gunzip_pet, [("in_files",         "in_file")]),
                (gunzip_pet,  warp,       [("out_file",         "image_to_align")]),

                # output
                (warp,        pet_output, [("normalized_image", "warped_files")]),
               ])

    return wf


def attach_spm_pet_preprocessing(main_wf, wf_name="spm_pet_preproc"):
    """ Attach a FDG-PET only pre-processing workflow that uses SPM12 to `main_wf`.

    Parameters
    ----------
    main_wf: nipype Workflow

    wf_name: str
        Name of the preprocessing workflow

    Nipype Inputs for `main_wf`
    ---------------------------
    Note: The `main_wf` workflow is expected to have an `input_files` and a `datasink` nodes.

    input_files.select.pet: input node

    datasink: nipype Node

    Returns
    -------
    main_wf: nipype Workflow
    """
    # Dependency workflows
    in_files = get_input_node(main_wf)
    datasink = get_datasink  (main_wf)

    # The base name of the 'pet' file for the substitutions
    pet_fbasename = remove_ext(op.basename(get_input_file_name(in_files, 'pet')))

    # get the PET preprocessing pipeline
    pet_wf = spm_pet_preproc(wf_name=wf_name)

    # dataSink output substitutions
    regexp_subst = [
                     (r"/w{pet}.nii", "/{pet}_mni.nii"),
                   ]
    regexp_subst = format_pair_list(regexp_subst, pet=pet_fbasename)
    regexp_subst += extension_duplicates(regexp_subst)
    datasink.inputs.regexp_substitutions = extend_trait_list(datasink.inputs.regexp_substitutions,
                                                             regexp_subst)

    # Connect the nodes
    main_wf.connect([
                # pet file input
                (in_files, pet_wf, [("pet", "pet_input.in_files")]),

                (pet_wf, datasink, [
                                    ("pet_output.warped_files",  "pet.@warped"),
                                   ]),
              ])

    return main_wf
