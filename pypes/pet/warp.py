# -*- coding: utf-8 -*-
"""
PET-only image registration nipype workflow.
"""
import os.path as op

from   ..preproc import spm_warp_to_mni
from   .._utils  import format_pair_list
from   ..utils   import (get_datasink,
                         extend_trait_list,
                         get_input_node,
                         remove_ext,
                         get_input_file_name,
                         extension_duplicates)


def attach_spm_pet_preprocessing(main_wf, wf_name='spm_pet_preproc'):
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
    warp_pet_wf = spm_warp_to_mni(wf_name=wf_name)

    # dataSink output substitutions
    regexp_subst = [
                    (r'/w{pet}.nii', '/{pet}_mni.nii'),
                   ]
    regexp_subst = format_pair_list(regexp_subst, pet=pet_fbasename)
    regexp_subst += extension_duplicates(regexp_subst)
    datasink.inputs.regexp_substitutions = extend_trait_list(datasink.inputs.regexp_substitutions,
                                                             regexp_subst)

    # Connect the nodes
    main_wf.connect([
                     # pet file input
                    (in_files,    warp_pet_wf, [('pet',                      'warp_input.in_files')]),
                    (warp_pet_wf, datasink,    [('warp_output.warped_files', 'pet.@warped'),]),
                   ])

    return main_wf
