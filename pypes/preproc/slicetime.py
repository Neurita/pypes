# -*- coding: utf-8 -*-
"""
Helper functions for Slice Timing Correction
"""
import os.path as op

import nipype.pipeline.engine as pe
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
    import nipype.interfaces.afni as afni

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
    in_files: str or list of str
        Path to the input file(s).

    out_prefix: str
        Prefix to the output file.
        Default: 'a'

    num_slices: int

    time_repetition: int or str
        The time repetition (TR) of the input dataset in seconds
        Default: 0
        If left to default will read the TR from the nifti image header.

    time_acquisition: int
        Time of volume acquisition. usually calculated as TR-(TR/num_slices)

    ignore_first: int
        Number of acquisitions thrown away from the beginning of the
        dataset.

    ref_slice: int
        Index of the reference slice

    slice_order: list of int

    Returns
    -------
    stc: nipype interface
    """
    import nipype.interfaces.spm  as spm

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


def nipy_fmrirealign4d(in_files=traits.Undefined,
                       time_repetition=2,
                       slice_order=None,
                       loops=5,
                       ):
    """ Return a nipype interface to the nipy.FmriRealign4d

    Parameters
    ----------
    in_files: str or list of str
        Path to the input file(s).

    time_repetition: int or str
        The time repetition (TR) of the input dataset in seconds
        Default: 0
        If left to default will read the TR from the nifti image header.

    slice_order: list of int
        This would be equivalent to
        entering np.argsort(spm_slice_order) for this field.
        This field will be deprecated in future Nipy releases and be
        replaced by actual slice acquisition
        times.

    loops: int
        Number of loops used to realignment runs.

    Returns
    -------
    stc: nipype interface
    """
    import nipype.interfaces.nipy as nipy

    stc = nipy.FmriRealign4d()

    stc.inputs.in_files         = in_files
    stc.inputs.time_repetition  = time_repetition
    stc.inputs.slice_order      = slice_order
    stc.inputs.loops            = loops

    #res = stc.run()
    return stc


def auto_spm_slicetime(in_files=traits.Undefined,
                       out_prefix=traits.Undefined,
                       num_slices=0,
                       time_repetition=-1,
                       time_acquisition=-1,
                       ignore_first=6,
                       ref_slice=traits.Undefined,
                       slice_order=None,
                       wf_name='auto_spm_slicetime'):
    """ A workflow that tries to automatically read the slice timing correction parameters
    from the input files and passes them to a spm.SliceTiming node.

    Parameters
    ----------
    in_files: str or list of str
        Path to the input file(s).

    out_prefix: str
        Prefix to the output file.
        Default: 'a'

    num_slices: int
        Number of slices of `in_files`.

    time_repetition: int or str
        The time repetition (TR) of the input dataset in seconds
        Default: 0
        If left to default will read the TR from the nifti image header.

    time_acquisition: int
        Time of volume acquisition. usually calculated as TR-(TR/num_slices)

    ignore_first: int
        Number of acquisitions thrown away from the beginning of the
        dataset.

    ref_slice: int
        Index of the reference slice

    slice_order: list of int
        List of integers with the order in which slices are acquired

    wf_name: str
        Name of the workflow

    Nipype Inputs
    -------------
    ## Mandatory:
    params.in_files:

    ## Optional:
    params.num_slices

    params.slice_order

    params.time_repetition

    params.time_acquisition

    params.ref_slice

    params.slice_mode

    Nipype Outputs
    --------------
    slice_timer.timecorrected_files

    Returns
    -------
    auto_spm_stc: nipype Workflow
        SPM slice timing correction workflow with automatic
        parameters detection.
    """
    def sum_one(slice_order): # SPM starts count from 1
        return [i+1 for i in slice_order]

    # Declare the processing nodes
    params = pe.Node(slice_timing_params(), name='params')
    stc    = pe.Node(spm_slicetime(in_files=in_files,
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
                (params, stc, [("in_files",                 "in_files"),
                               ("num_slices",               "num_slices"),
                               ("ref_slice",                "ref_slice"),
                               (("slice_order", sum_one),   "slice_order"),
                               ("time_acquisition",         "time_acquisition"),
                               ("time_repetition",          "time_repetition"),
                              ]),
              ])

    return wf


def auto_nipy_slicetime(in_files=traits.Undefined,
                        time_repetition=-1,
                        slice_order=None,
                        loops=5,
                        wf_name='auto_spm_slicetime'):
    """ A workflow that tries to automatically read the slice timing correction parameters
    from the input files and passes them to a nipy.fMRIRealign4D node.

    Parameters
    ----------
    in_files: str or list of str
        Path to the input file(s).

    time_repetition: int or str
        The time repetition (TR) of the input dataset in seconds
        Default: 0
        If left to default will read the TR from the nifti image header.

    slice_order: list of int
        List of integers with the order in which slices are acquired

    loops: int
        Number of loops used to realignment runs.

    wf_name: str
        Name of the workflow

    Nipype Inputs
    -------------
    ## Mandatory:
    params.in_files:

    params.time_repetition

    ## Optional:
    params.slice_order

    params.ref_slice

    params.slice_mode

    Nipype Outputs
    --------------
    slice_timer.out_file

    slice_timer.par_file

    Returns
    -------
    auto_nipy_stc: nipype Workflow
        Nipy 4D alignment and slice timing correction workflow with automatic
        parameters detection.
    """
    # Declare the processing nodes
    params = pe.Node(slice_timing_params(), name='params')
    stc    = pe.Node(nipy_fmrirealign4d(in_files=in_files,
                                        time_repetition=time_repetition,
                                        slice_order=slice_order,
                                        loops=loops),
                     name='slice_timer')

    # Create the workflow object
    wf = pe.Workflow(name=wf_name)

    # Connect the nodes
    wf.connect([
                (params, stc, [("in_files",         "in_files"),
                               ("slice_order",      "slice_order"),
                               ("time_repetition",  "time_repetition"),
                              ]),
              ])

    return wf