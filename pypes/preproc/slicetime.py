"""
Helper functions for Slice Timing Correction
"""
import os.path as op

import nibabel as nib
import nipype.interfaces.afni as afni
import nipype.interfaces.spm as spm
from   nipype.interfaces.base import traits, isdefined

from ..utils import remove_ext
from .._utils import check_equal


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


def spm_siemens_slicetime(in_file=traits.Undefined,
                          out_prefix=traits.Undefined,
                          tr=-1,
                          time_acquisition=-1,
                          ignore_first=traits.Undefined,
                          ref_slice=traits.Undefined,
                          slice_mode='unknown',
                          ):
    """ Return a nipype interface to the SiemensSliceTiming correction which
    helps setting some parameters for the SPM SliceTiming correction specially
    regarding Siemens acquisitions.

    Parameters
    ----------
    in_file: str
        Path to the input file

    out_prefix: str
        Prefix to the output file.
        Default: 'a'

    tr: int or str
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

    slice_mode: str
        Choices:
            'unknown': auto detect if images are from Siemens and converted with dcm2nii from Nov 2013 or later #kNIFTI_SLICE_UNKNOWN
            'seq_inc': sequential ascending kNIFTI_SLICE_SEQ_INC = 1; %1,2,3,4
            'seq_dec': sequential descending kNIFTI_SLICE_SEQ_DEC = 2; %4,3,2,1
            'alt_inc': Siemens: interleaved ascending with odd number of slices, interleaved for other vendors kNIFTI_SLICE_ALT_INC = 3; %1,3,2,4
            'alt_dec': descending interleaved kNIFTI_SLICE_ALT_DEC = 4; %4,2,3,1
            'alt_inc2': Siemens interleaved ascending with even number of slices kNIFTI_SLICE_ALT_INC2 = 5; %2,4,1,3
            'alt_dec2': Siemens interleaved descending with even number of slices kNIFTI_SLICE_ALT_DEC2 = 6; %3,1,4,2

        Default: 'unknown'
        If left to default will read the TR from the nifti image header.

    Returns
    -------
    stc: nipype interface
    """
    stc = SiemensSliceTiming()

    stc.slice_mode = slice_mode
    stc.inputs.in_file = in_file
    stc.inputs.out_prefix = out_prefix
    stc.inputs.out_prefix = out_prefix
    stc.inputs.ref_slice = ref_slice
    stc.inputs.ignore = ignore_first

    # the automatically calculated settings
    stc.inputs.slice_order = [0]
    stc.inputs.time_repetition = -1
    stc.inputs.time_acquisition = -1
    stc.inputs.ref_slice = 0
    stc.inputs.num_slices = 0

    #res = stc.run()
    return stc