# -*- coding: utf-8 -*-
"""
PET-only image registration nipype workflow.
"""
import os.path as op

import nipype.pipeline.engine    as pe
import nipype.interfaces.spm     as spm
from   nipype.interfaces.utility import IdentityInterface
from   nipype.algorithms.misc    import Gunzip

from   ..preproc import spm_normalize
from   ..config  import setup_node
from   .._utils  import format_pair_list
from   ..utils   import (get_datasink,
                         extend_trait_list,
                         get_input_node,
                         remove_ext,
                         get_input_file_name,
                         extension_duplicates)


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
                                             name="pet_input",)

    gunzip_pet  = setup_node(Gunzip(), name="gunzip_pet",)

    warp = setup_node(spm.Normalize12(jobtype='estwrite',
                                      affine_regularization_type='mni'),
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
