# -*- coding: utf-8 -*-
"""
Nipype workflows to process diffusion MRI.
"""
import os.path as op

import nipype.interfaces.spm     as spm
import nipype.pipeline.engine    as pe
from   nipype.interfaces.fsl     import ExtractROI, BET, EddyCorrect, DTIFit
from   nipype.interfaces.io      import DataSink, SelectFiles

from .utils import find_wf_node


def fsl_dti_preprocessing(wf_name="fsl_dti_preproc"):
    """ Run the diffusion MRI pre-processing workflow against the diff files in `data_dir`.

    It does:
    - EddyCorrect
    - DTIFit

    Nipype Inputs
    -------------
    eddy_correct.in_file: traits.File
        path to the diffusion MRI image
    roi.in_file: traits.File
        path to the diffusion MRI image
    dtifit.bvecs: traits.File
        path to the b vectors file
    dtifit.bvals: traits.File
        path to the b values file

    Returns
    -------
    wf: nipype Workflow
    """
    roi          = pe.Node(ExtractROI(),            name="roi")
    bet          = pe.Node(BET(),                   name="bet")
    eddy_correct = pe.Node(EddyCorrect(),           name="eddy_correct")
    dtifit       = pe.Node(DTIFit(),                name="dtifit")

    roi.inputs.t_min = 0
    roi.inputs.t_size = 1

    bet.inputs.mask = True

    # Create the workflow object
    wf = pe.Workflow(name=wf_name)

    # Connect the nodes
    wf.connect([
                (roi,          bet,         [("roi_file", "in_file")]),
                (eddy_correct, dtifit,      [("eddy_corrected", "dwi")]),
                (bet,          dtifit,      [("mask_file", "mask")]),
              ])
    return wf


def attach_fsl_dti_preprocessing(main_wf, wf_name="fsl_dti_preproc"):
    """ Attach the FSL-based diffusion MRI pre-processing workflow to the `main_wf`.

    Parameters
    ----------
    main_wf: nipype Workflow

    wf_name: str
        Name of the preprocessing workflow

    Nipype Inputs for `main_wf`
    ---------------------------
    Note: The `main_wf` workflow is expected to have an `input_files` and a `datasink` nodes.

    input_files.select.diff: input node

    datasink: nipype Node

    Returns
    -------
    main_wf: nipype Workflow
    """
    in_files = find_wf_node(main_wf, SelectFiles)
    datasink = find_wf_node(main_wf, DataSink)

    # The workflow box
    dti_wf = fsl_dti_preprocessing(wf_name=wf_name)

    # input and output diffusion MRI workflow to main workflow connections
    main_wf.connect([(in_files, dti_wf,   [("diff",                                  "eddy_correct.in_file"),
                                           ("diff",                                  "roi.in_file"),
                                           ("diff_bval",                             "dtifit.bvals"),
                                           ("diff_bvec",                             "dtifit.bvecs")]),
                     (dti_wf,   datasink, [("eddy_correct.eddy_corrected",           "diff.@eddy_corrected"),
                                           ("dtifit.V1",                             "diff.@v1"),
                                           ("dtifit.V2",                             "diff.@v2"),
                                           ("dtifit.V3",                             "diff.@v3"),
                                           ("dtifit.L1",                             "diff.@l1"),
                                           ("dtifit.L2",                             "diff.@l2"),
                                           ("dtifit.L3",                             "diff.@l3"),
                                           ("dtifit.MD",                             "diff.@mean_diffusivity"),
                                           ("dtifit.FA",                             "diff.@fractional_anisotropy"),
                                           ("dtifit.MO",                             "diff.@mode_of_anisotropy"),
                                           ("dtifit.S0",                             "diff.@s0"),
                                           ("dtifit.tensor",                         "diff.@tensor"),],),
                    ])

    return main_wf
