"""
Helper functions for Slice Timing Correction
"""
import os.path as op

from ..utils import remove_ext

import nipype.interfaces.afni as afni
from   nipype.interfaces.base import traits, isdefined


def afni_slicetime(in_file=traits.Undefined,
                   out_file=traits.Undefined,
                   tr=traits.Undefined,
                   tr_units='s',
                   ignore_first=traits.Undefined,
                   slice_mode=traits.Undefined,
                   tslice=traits.Undefined,
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

    tslice: int
        Align each slice to the time offset of the given `tslice` index.
        Note: only one of the tzero and slice options can be used.

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
        tshift.inputs.tr = '{}{}'.format(tr, tr_units)

    if isdefined(in_file) and not isdefined(out_file):
        out_file = op.join(op.dirname(in_file), 'tshift_' + remove_ext(op.basename(in_file)))

    tshift.inputs.in_file = in_file
    tshift.inputs.out_file = out_file
    tshift.inputs.tpattern = slice_mode
    tshift.inputs.tzero = tzero
    tshift.inputs.tslice = tslice
    tshift.inputs.ignore = ignore_first
    tshift.inputs.outputtype = out_type

    #res = tshift.run()
    return tshift
