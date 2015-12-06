"""
Helper functions for Slice Timing Correction
"""
import os.path as op

import nipype.pipeline.engine as pe
import nipype.interfaces.afni as afni
import nipype.interfaces.spm  as spm
from   nipype.interfaces.base import traits, isdefined

from   .slicetime_params import slice_timing_params
from   ..utils import remove_ext


def afni_slicetime(in_file=traits.Undefined,
                   out_file=traits.Undefined,
                   tr=traits.Undefined,
                   tr_units='s',
                   ignore_first=traits.Undefined,
                   slice_mode=traits.Undefined,
                   ref_slice=traits.Undefined,
                   tzero=traits.Undefined,
                   out_type='NIFTI_GZ'):
    """ Return a nipype interface to the AFNI 3dTshift command.

    Parameters
    ----------
    in_file: str
        Path to the input file

    out_file: str
        Path to the output file.

    tr: int or str
        The TR of the input dataset

    tr_units: str
        The units of `tr`. Choices: ('s', 'ms')

    ignore_first: int
        Number of acquisitions thrown away from the beginning of the
        dataset.

    tzero: float
        Align each slice to time offset `tzero`.
        The value of 'tzero' must be between the
        minimum and maximum slice temporal offsets.
        Note: The default alignment time is the average
              of the 'tpattern' values (either from the
              dataset header or from the -tpattern option)

    ref_slice: int
        Align each slice to the time offset of the given `ref_slice` index.
        Note: only one of the tzero and slice options can be used.
        -tslice flag

    slice_mode: str
        Choices:
        'alt+z' or 'altplus'    = alternating in the plus direction
        'alt+z2'                = alternating, starting at slice #1 instead of #0
        'alt-z' or 'altminus'   = alternating in the minus direction
        'alt-z2'                = alternating, starting at slice #nz-2 instead of #nz-1
        'seq+z' or 'seqplus'    = sequential in the plus direction
        'seq-z' or 'seqminus'   = sequential in the minus direction
         @filename              = read temporal offsets from 'filename'

    out_type: str
        ('NIFTI_GZ' or 'AFNI' or 'NIFTI')
        AFNI output filetype

    Returns
    -------
    tshift: nipype interface
    """
    tshift = afni.TShift()

    if isdefined(tr):
        tr = '{}{}'.format(tr, tr_units)

    if isdefined(in_file):
        if not isdefined(out_file):
            out_file = op.join(op.dirname(in_file), 'tshift_' + remove_ext(op.basename(in_file)))

    tshift.inputs.in_file = in_file
    tshift.inputs.out_file = out_file
    tshift.inputs.tpattern = slice_mode
    tshift.inputs.tzero = tzero
    tshift.inputs.tr = tr
    tshift.inputs.tslice = ref_slice
    tshift.inputs.ignore = ignore_first
    tshift.inputs.outputtype = out_type

    #res = tshift.run()
    return tshift


def spm_slicetime(in_files=traits.Undefined,
                  out_prefix=traits.Undefined,
                  num_slices=0,
                  time_repetition=-1,
                  time_acquisition=-1,
                  ignore_first=traits.Undefined,
                  ref_slice=traits.Undefined,
                  slice_order=None,
                  ):
    """ Return a nipype interface to the SPM SliceTiming.

    Parameters
    ----------
    in_files: str
        Path to the input file

    out_prefix: str
        Prefix to the output file.
        Default: 'a'

    num_slices: int

    time_repetition: int or str
        The time repetition (TR) of the input dataset in seconds
        Default: 0
        If left to default will read the TR from the nifti image header.

    ignore_first: int
        Number of acquisitions thrown away from the beginning of the
        dataset.

    tzero: float
        Align each slice to time offset `tzero`.
        The value of 'tzero' must be between the
        minimum and maximum slice temporal offsets.
        Note: The default alignment time is the average
              of the 'tpattern' values (either from the
              dataset header or from the -tpattern option)

    ref_slice: int
        Index of the reference slice

    slice_order: list of int

    Returns
    -------
    stc: nipype interface
    """
    stc = spm.SliceTiming()

    stc.inputs.in_files         = in_files
    stc.inputs.out_prefix       = out_prefix
    stc.inputs.slice_order      = slice_order
    stc.inputs.ref_slice        = ref_slice
    stc.inputs.ignore           = ignore_first
    stc.inputs.time_repetition  = time_repetition
    stc.inputs.time_acquisition = time_acquisition
    stc.inputs.num_slices       = num_slices

    #res = stc.run()
    return stc


def auto_spm_slicetime(in_files=traits.Undefined,
                       out_prefix=traits.Undefined,
                       num_slices=0,
                       time_repetition=-1,
                       time_acquisition=-1,
                       ignore_first=traits.Undefined,
                       ref_slice=traits.Undefined,
                       slice_order=None,
                       wf_name='auto_spm_slicetime'):
    """

    Parameters
    ----------
    in_files
    out_prefix
    num_slices
    time_repetition
    time_acquisition
    ignore_first
    ref_slice
    slice_order
    wf_name

    Nipype Inputs
    -------------
    - Mandatory:
    params.in_files:


    - Optional:
    params.num_slices

    params.slice_order

    params.time_repetition

    params.time_acquisition

    params.ref_slice

    params.slice_mode

    Nipype Outputs
    --------------
    slice_timer.out_files

    Returns
    -------
    auto_spm_stc: nipype Workflow
    SPM slice timing correction workflow with automatic
    parameters detection
    """

    params = pe.Node(slice_timing_params(), name='params')
    stc = pe.Node(spm_slicetime(in_files=in_files,
                                out_prefix=out_prefix,
                                num_slices=num_slices,
                                time_repetition=time_repetition,
                                time_acquisition=time_acquisition,
                                ignore_first=ignore_first,
                                ref_slice=ref_slice,
                                slice_order=slice_order), name='slice_timer')

        # Create the workflow object
    wf = pe.Workflow(name=wf_name)

    # Connect the nodes
    wf.connect([
                (params, stc, [("in_files",         "in_files"),
                               ("num_slices",       "num_slices"),
                               ("slice_order",      "slice_order"),
                               ("time_repetition",  "time_repetition"),
                               ("time_acquisition", "time_acquisition"),
                               ("ref_slice",        "ref_slice"),
                              ]),
              ])

    return wf