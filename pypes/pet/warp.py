# -*- coding: utf-8 -*-
"""
PET-only image registration nipype workflow.
"""
import os.path as op

import nipype.pipeline.engine    as pe
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
                       get_input_file_name)

from   .._utils import (flatten_list,
                        format_pair_list)

def printit(x):
    print(x)
    return x


def spm_pet_preprocessing(wf_name="spm_pet_preproc"):
    """ Run a PET-only pre-processing workflow against the gunzip_pet.in_file files.
    It depends on the anat_preproc_workflow, so if this has not been run, this function
    will run it too.

    It does:
    - Warp each individual PET image to the default (SPM) PET template (H2O),
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
    pet_input.in_files: list of traits.File
        The raw NIFTI_GZ PET image files

    pet_input.tissues: list of traits.File

    Nipype outputs
    --------------
    pet_output.warped_files: list of existing file
        The warped PET files.

    pet_output.template: existing file
        The common custom PET template file.

    pet_output.corrected_files: list of existing file
        The warped PVC PET files

    Returns
    -------
    wf: nipype Workflow
    """
    # input
    pet_input = pe.MapNode(IdentityInterface(fields=["in_files"]),
                                             name="pet_input",
                                             iterfield=["in_files"]
                                             )

    # coreg pet
    gunzip_pet  = setup_node(Gunzip(),
                             name="gunzip_pet",
                             type='map',
                             iterfield=['in_file'])

    warp2mni    = setup_node(spm_normalize(voxel_size=[2, 2, 2]),
                             type='map',
                             iterfield=['image_to_align'],
                             name="warp_pet2mni")

    pet_warped = pe.JoinNode(interface=IdentityInterface(fields=["in_files"]),
                             joinsource="warp_pet2mni",
                             joinfield="in_files", name="pet_warped")

    # output
    pet_output = setup_node(IdentityInterface(fields=["warped_files",
                                                     ]),
                                               name="pet_output")

    print_out = pe.Node(Function(input_names=['x'], output_names=['x'], function=printit),
                        name='print')


    # Create the workflow object
    wf = pe.Workflow(name=wf_name)

    wf.connect([
                # inputs
                (pet_input,   gunzip_pet, [("in_files",         "in_file")]),
                (gunzip_pet,  warp2mni,   [("out_file",         "image_to_align")]),

                (warp2mni,  pet_warped,   [("normalized_image", "in_files")]),

                (pet_warped,  print_out,   [("in_files", "x")]),
                # output
                (warp2mni,    pet_output, [("normalized_image", "warped_files")]),
               ])

    return wf


def attach_spm_pet_preprocessing(main_wf, wf_name="spm_pet_preproc"):
    """ Attach a FDG-PET only pre-processing workflow that uses SPM12 to `main_wf`.

    Nipype Inputs for `main_wf`
    ---------------------------
    Note: The `main_wf` workflow is expected to have an `input_files` and a `datasink` nodes.

    input_files.select.pet: input node

    datasink: nipype Node

    Parameters
    ----------
    main_wf: nipype Workflow

    wf_name: str
        Name of the preprocessing workflow

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

    # The base name of the 'pet' file for the substitutions
    pet_fbasename = remove_ext(op.basename(get_input_file_name(in_files, 'pet')))

    # get the PET preprocessing pipeline
    pet_wf = spm_pet_preprocessing(wf_name=wf_name)

    # dataSink output substitutions
    # regexp_subst = [
    #                ]
    # regexp_subst = format_pair_list(regexp_subst, pet=pet_fbasename)
    # datasink.inputs.regexp_substitutions = extend_trait_list(datasink.inputs.regexp_substitutions,
    #                                                          regexp_subst)

    # Connect the nodes
    main_wf.connect([
                # pet file input
                (in_files, pet_wf, [("pet", "pet_input.in_files")]),

                (pet_wf, datasink, [
                                    ("pet_output.warped_files",  "pet.@warped"),
                                   ]),
              ])

    return main_wf
