# -*- coding: utf-8 -*-
"""
Nipype workflows to process anatomical MRI.
"""
import os.path as op

import nipype.interfaces.spm     as spm
import nipype.pipeline.engine    as pe
from   nipype.algorithms.misc    import Gunzip, TSNR
from   nipype.interfaces.utility import Function, Select, Split, Merge, IdentityInterface
from   nipype.interfaces.nipy.preprocess import Trim
from   nipype.interfaces import fsl

from   ..preproc import (spm_apply_deformations,
                         auto_nipy_slicetime,
                         auto_spm_slicetime,
                         nipy_motion_correction,
                         extract_noise_components,
                         )

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
    rest_input.in_files: traits.File
        path to the resting-state image

    Nipype Outputs
    --------------
    rest_output.out_file: traits.File

    rest_output.motion_corrected: traits.File
        The motion corrected file.

    rest_output.motion_params: traits.File
        The affine transformation file.

    Returns
    -------
    wf: nipype Workflow
    """
    # Create the workflow object
    wf = pe.Workflow(name=wf_name)

    # input identities
    rest_input = setup_node(IdentityInterface(fields=["in_files",
                                                      "num_noise_components",
                                                      "highpass_sigma",
                                                      "lowpass_sigma",
                                                      ], mandatory_inputs=True),
                            name="rest_input")

    # rs-fMRI preprocessing nodes
    trim    = setup_node(Trim(),
                         name="trim")
    stc_wf  = auto_spm_slicetime()
    realign = nipy_motion_correction()



    # output identities
    rest_output = setup_node(IdentityInterface(fields=["out_file",
                                                       "motion_corrected",
                                                       "motion_params",
                                                       ],
                                               mandatory_inputs=True),
                             name="rest_output")


    # Connect the nodes
    wf.connect([
                # trim
                (rest_input,   trim,    [("in_files",  "in_file")]),

                #slice time correction
                (trim,        stc_wf,   [("out_file",  "stc_params.in_files")]),

                # motion correction
                (stc_wf,      realign,  [("timecorrected_files", "in_file")]),

                # motion statistics

                # nuisance correction

                # median angle?

                # temporal filtering

                # registration

                #output test
                (trim,   rest_output,  [("out_file",  "out_file")]),
                (stc_wf, rest_output,  [("out_file",  "motion_corrected")]),
                (stc_wf, rest_output,  [("par_file",  "motions_params")]),
                # output
              ])
    return wf


def trash():
    inputnode = pe.Node(interface=util.IdentityInterface(fields=['func',
                                                                 'num_noise_components',
                                                                 'highpass_sigma',
                                                                 'lowpass_sigma'
                                                                 ]),
                        name='inputspec')

    outputnode = pe.Node(interface=util.IdentityInterface(fields=[
        'noise_mask_file',
        'filtered_file',
    ]),
        name='outputspec')
    slicetimer = pe.Node(fsl.SliceTimer(), name='slicetimer')
    realigner = create_realign_flow()
    tsnr = pe.Node(TSNR(regress_poly=2), name='tsnr')
    getthresh = pe.Node(interface=fsl.ImageStats(op_string='-p 98'),
                        name='getthreshold')
    threshold_stddev = pe.Node(fsl.Threshold(), name='threshold')
    compcor = pe.Node(util.Function(input_names=['realigned_file',
                                                 'noise_mask_file',
                                                 'num_components'],
                                    output_names=['noise_components'],
                                    function=extract_noise_components),
                      name='compcorr')
    remove_noise = pe.Node(fsl.FilterRegressor(filter_all=True),
                           name='remove_noise')
    bandpass_filter = pe.Node(fsl.TemporalFilter(),
                              name='bandpass_filter')

    # Define connections
    restpreproc.connect(inputnode, 'func', slicetimer, 'in_file')
    restpreproc.connect(slicetimer, 'slice_time_corrected_file',
                        realigner, 'inputspec.func')
    restpreproc.connect(realigner, 'outputspec.realigned_file', tsnr, 'in_file')
    restpreproc.connect(tsnr, 'stddev_file', threshold_stddev, 'in_file')
    restpreproc.connect(tsnr, 'stddev_file', getthresh, 'in_file')
    restpreproc.connect(getthresh, 'out_stat', threshold_stddev, 'thresh')
    restpreproc.connect(realigner, 'outputspec.realigned_file',
                        compcor, 'realigned_file')
    restpreproc.connect(threshold_stddev, 'out_file',
                        compcor, 'noise_mask_file')
    restpreproc.connect(inputnode, 'num_noise_components',
                        compcor, 'num_components')
    restpreproc.connect(tsnr, 'detrended_file',
                        remove_noise, 'in_file')
    restpreproc.connect(compcor, 'noise_components',
                        remove_noise, 'design_file')
    restpreproc.connect(inputnode, 'highpass_sigma',
                        bandpass_filter, 'highpass_sigma')
    restpreproc.connect(inputnode, 'lowpass_sigma',
                        bandpass_filter, 'lowpass_sigma')
    restpreproc.connect(remove_noise, 'out_file', bandpass_filter, 'in_file')
    restpreproc.connect(threshold_stddev, 'out_file',
                        outputnode, 'noise_mask_file')
    restpreproc.connect(bandpass_filter, 'out_file',
                        outputnode, 'filtered_file')
    return restpreproc


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
    main_wf.connect([(in_files, rest_wf,  [("rest",                   "rest_input.in_files")]),

                     # test output
                     (rest_wf,  datasink, [("rest_output.out_file",   "rest.@trim")],),
                    ])

    return main_wf
