# -*- coding: utf-8 -*-
"""
Nipype workflows to process anatomical MRI.
"""
import os.path as op

import nipype.interfaces.spm     as spm
import nipype.pipeline.engine    as pe
from   nipype.algorithms.misc    import Gunzip
from   nipype.interfaces.utility import Function, Select, Split, Merge, IdentityInterface
from nipype.interfaces.nipy.preprocess import Trim

from   ..preproc import spm_apply_deformations, slice_timing_params

from   .._utils import format_pair_list
from   ..utils import (setup_node,
                       remove_ext,
                       extend_trait_list,
                       get_input_node,
                       get_datasink,
                       get_input_file_name)



def rest_preprocessing_wf(wf_name="rest_preproc"):
    """ Run the resting-state fMRI pre-processing workflow against the rest files in `data_dir`.

    It does:
    - Trim first 6 volumes of the rs-fMRI file.
    - Slice Timing Correction
    -

    Nipype Inputs
    -------------
    rest_input.in_file: traits.File
        path to the resting-state image


    Nipype Outputs
    --------------
    rest_output.out_file: traits.File


    Returns
    -------
    wf: nipype Workflow
    """
    # input identities
    rest_input = pe.Node(IdentityInterface(fields=["in_file"], mandatory_inputs=True), #, "tissues", "anat", "mni_to_anat"],
                         name="rest_input")

    # rs-fMRI preprocessing nodes
    trim      = setup_node(Trim(), name="trim")
    st_params = pe.Node(slice_timing_params(), name='st_params')

    # output identities
    rest_output = pe.Node(IdentityInterface(fields=["out_file"], mandatory_inputs=True),
                          name="rest_output")

    # Create the workflow object
    wf = pe.Workflow(name=wf_name)

    # Connect the nodes
    wf.connect([
                # trim
                (rest_input,    trim,    [("rest",     "in_file")]),

                #slice time correction
                (trim,     st_params,   [("out_file",  "in_file")]),

                #output test
                (trim, rest_output,     [("out_file",  "out_file")]),

                # output
              ])
    return wf


def attach_rest_preprocessing(main_wf, wf_name="rest_preproc", params=None):
    """ Attach the resting-state MRI pre-processing workflow to the `main_wf`.

    Parameters
    ----------
    main_wf: nipype Workflow

    wf_name: str
        Name of the preprocessing workflow

    Nipype Inputs for `main_wf`
    ---------------------------
    Note: The `main_wf` workflow is expected to have an `input_files` and a `datasink` nodes.

    input_files.select.anat: input node

    datasink: nipype Node

    Returns
    -------
    main_wf: nipype Workflow
    """
    in_files = get_input_node(main_wf)
    datasink = get_datasink  (main_wf)

    # The workflow box
    rest_wf = rest_preprocessing_wf(wf_name=wf_name)

    # The base name of the 'rest' file for the substitutions
    rest_fbasename = remove_ext(op.basename(get_input_file_name(in_files, 'rest')))

    # dataSink output substitutions
    regexp_subst = [
                   ]
    regexp_subst = format_pair_list(regexp_subst, anat=rest_fbasename)
    datasink.inputs.regexp_substitutions = extend_trait_list(datasink.inputs.regexp_substitutions,
                                                             regexp_subst)

    # input and output anat workflow to main workflow connections
    main_wf.connect([(in_files, rest_wf,  [("rest",                   "rest_input.in_file")]),

                     # test output
                     (rest_wf,  datasink, [("rest_output.out_file",   "rest.@trim")],),
                    ])

    return main_wf
