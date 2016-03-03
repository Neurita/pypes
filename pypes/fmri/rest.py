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

from   pypes.preproc import spm_apply_deformations
from   pypes._utils import format_pair_list
from   pypes.utils import (remove_ext,
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
    rest_input.rest: traits.File
        path to the resting-state image

    rest_input.tissues: traits.File
        path to the

    Returns
    -------
    wf: nipype Workflow
    """
    dti_input    = pe.Node(IdentityInterface(fields=["rest"], mandatory_inputs=True), #, "tissues", "anat", "mni_to_anat"],
                           name="rest_input")

    # rs-fMRI preprocessing nodes
    trim        = pe.Node(Trim(begin_idx=6),        name="trim")
    slice_time  =


    # Create the workflow object
    wf = pe.Workflow(name=wf_name)

    # Connect the nodes
    wf.connect([
                # new segment
                (biascor,      gunzip_anat, [("output_image", "in_file"      )]),
                (gunzip_anat,  segment,     [("out_file",     "channel_files")]),

                # Normalize12
                (segment, warp_anat, [("forward_deformation_field", "deformation_file")]),
                (segment, warp_anat, [("bias_corrected_images",     "apply_to_files")]),
              ])
    return wf


def attach_rest_preprocessing(main_wf, wf_name="rest_preproc"):
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
    rest_wf = rest_preprocessing(wf_name=wf_name)

    # The base name of the 'rest' file for the substitutions
    rest_fbasename = remove_ext(op.basename(get_input_file_name(in_files, 'rest')))

    # dataSink output substitutions
    regexp_subst = [
                   ]
    regexp_subst = format_pair_list(regexp_subst, anat=rest_fbasename)
    datasink.inputs.regexp_substitutions = extend_trait_list(datasink.inputs.regexp_substitutions,
                                                             regexp_subst)

    # input and output anat workflow to main workflow connections
    main_wf.connect([(in_files, rest_wf,  [("rest",                                  "input_image")]),
                     (rest_wf,  datasink, [("warp_anat.normalized_files",            "anat.@mni")],),
                     (rest_wf,  datasink, [("new_segment.modulated_class_images",    "anat.tissues.@warped"),
                                           ("new_segment.native_class_images",       "anat.tissues.@native"),
                                           ("new_segment.transformation_mat",        "anat.transform.@linear"),
                                           ("new_segment.forward_deformation_field", "anat.transform.@forward"),
                                           ("new_segment.inverse_deformation_field", "anat.transform.@inverse"),
                                           ("new_segment.bias_corrected_images",     "anat.@biascor"),
                                          ]),
                    ])

    return main_wf
