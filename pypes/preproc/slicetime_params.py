# -*- coding: utf-8 -*-
# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
"""
 Helper functions and nipype interface for
 reading slice timing correction parameters specially from Siemens acquisitions
"""
import os.path as op
import nibabel as nib
import numpy as np

# Nipype imports
from nipype.interfaces.utility import Function

from .._utils import check_equal, grep
from ..preproc.dicom import split_dcm_ahdr, dcm_ascii_hdr


def slicing_mode(dcm_file):
    """ Return the slicing mode of the fMRI acquisition file given
    one of its DICOM files. Avoid giving the first DICOM file.

    Parameters
    ----------
    dcm_file: str
       Path to the DICOM file

    Returns
    -------
    mode: str
        Choices: ('ascending', 'descending', 'interleaved')

    References
    ----------
    https://wiki.cimec.unitn.it/tiki-index.php?page=MRIBOLDfMRI
    http://www.mccauslandcenter.sc.edu/CRNL/tools/stc

    https://www.jiscmail.ac.uk/cgi-bin/webadmin?A2=ind03&L=SPM&D=0&P=5914721
    The timing can be found out most easily by looking into some shadow
    information: In any image from the syngo systems the measurement-protocol is
    present (in ascii). You may find it by looking for "### ASCCONV BEGIN ###"
    and "### ASCCONV END ###". Here you may also find the slice positions of the
    single images.
    > Now find "sSliceArray.ucMode". This is
    > 0x1 for ascending
    > 0x2 for descending
    > 0x4 for interleaved
    > In the interleaved mode, the information given by Peter Erhard is correct;
    for the rest it's clear anyway.

    Example
    -------
    import os.path as op
    from glob import glob
    dcmdir   = '/home/alexandre/data/pet_stadhauders/hiswork/FTD/'
    restdirs = [glob(op.join(op.abspath(sd), '*ep2d*')) for sd in glob(op.join(dcmdir, '*'))]
    dcms    = [glob(op.join(rd[0], '*.dcm'))[10] for rd in restdirs if rd]
    modes   = [slicing_mode(dcm) for dcm in dcms]
    print(modes)
    """
    _, ascconv = split_dcm_ahdr(dcm_ascii_hdr(dcm_file))

    mode_code = grep(ascconv, 'sSliceArray.ucMode')[0].split('=')[1].strip()

    code_modes = {'0x1': 'ascending',
                  '0x2': 'descending',
                  '0x4': 'interleaved'}

    return code_modes[mode_code]


class STCParameters(object):
    """ Class to calculate the parameters needed for slice timing correction.
    Some options are automated for Siemens acquisitions.

    Auto detection of slice order, i.e., slice_order == [0] and slice_mode == 'unknown'
    only works for images from Siemens and converted with dcm2nii from Nov 2013 or later.

    Siemens have unusual interleaving
    - http://cbs.fas.harvard.edu/node/559#slice_order
    - https://wiki.cimec.unitn.it/tiki-index.php?page=MRIBOLDfMRI

    This class is based on the script by Chris Rorden's Neuropsychology Lab in:
    - http://www.mccauslandcenter.sc.edu/CRNL/tools/stc

    This is a callable class so you can instance this class in a nipype Function object.
    See `slice_timing_params` in this module after the declaration of this class.

    Parameters
    ----------
    in_files: str
        Path to the input files

    num_slices: int
        Number of slices of `in_files`.

    ref_slice: int
        Index of the reference slice

    slice_order: list of ints
        List of integers with the order in which slices are acquired

    time_acquisition: int
        Time of volume acquisition. usually calculated as TR-(TR/num_slices)

    time_repetition: int
        The time repetition (TR) of the input dataset in seconds
        Default: 0
        If left to default will read the TR from the nifti image header.

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
        If left to default will try to detect the TR from the nifti image header, if it doesn't work
        an AttributeError exception will be raise.
    """
    def __init__(self, in_files=None,
                       num_slices=0,
                       time_repetition=None,
                       time_acquisition=None,
                       slice_order=None,
                       ref_slice=None,
                       slice_mode='unknown',):

        self.in_files         = in_files
        self.num_slices       = num_slices
        self.slice_order      = slice_order
        self.time_repetition  = time_repetition
        self.time_acquisition = time_acquisition
        self.ref_slice        = ref_slice
        self.slice_mode       = slice_mode

    def __call__(self, in_files,
                       num_slices=0,
                       ref_slice=None,
                       slice_order=None,
                       time_acquisition=None,
                       time_repetition=None,
                       slice_mode='unknown'):
        """

        Parameters
        ----------
        in_files: str
            Path to the input files

        num_slices: int
            Number of slices of `in_files`.

        ref_slice: int
            Index of the reference slice

        slice_order: list of ints
            List of integers with the order in which slices are acquired

        time_acquisition: int
            Time of volume acquisition. usually calculated as TR-(TR/num_slices)

        time_repetition: int
            The time repetition (TR) of the input dataset in seconds
            Default: 0
            If left to default will read the TR from the nifti image header.

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
            If left to default will try to detect the TR from the nifti image header, if it doesn't work
            an AttributeError exception will be raise.

        Returns
        -------
        num_slices

        ref_slice

        slice_order

        time_acquisition

        time_repetition
        """

        # If you have used this class for more than once and this exception is raised check the following comment:
        # TODO: decide to remove `isdefined` or change the `is not None` checks in the `set` functions
        from nipype.interfaces.base import isdefined
        if not isdefined(num_slices):
            raise ValueError("Expected `None` or an integer for `num_slices`, got a traits.Undefined.")

        self.in_files         = in_files
        self.num_slices       = num_slices
        self.time_repetition  = time_repetition
        self.time_acquisition = time_acquisition
        self.slice_order      = slice_order
        self.ref_slice        = ref_slice
        self.slice_mode       = slice_mode

        return self.fit()

    def fit(self):
        """

        Returns
        -------
        num_slices

        ref_slice

        slice_order

        time_acquisition

        time_repetition
        """
        num_slices       = self.set_num_slices()
        slice_order      = self.set_slice_order()
        time_repetition  = self.set_time_repetition()
        time_acquisition = self.set_time_acquisition()
        ref_slice        = self.set_ref_slice()

        return num_slices, ref_slice, slice_order, time_acquisition, time_repetition

    def _check_in_files(self):
        if isinstance(self.in_files, str):
            in_files = [self.in_files]
        else:
            in_files = self.in_files

        for f in in_files:
            if not op.exists(f):
                raise IOError('Expected an existing file in `in_files`')

        return in_files

    @staticmethod
    def _get_n_slices(in_file):
        img = nib.load(in_file)

        n_slices = 0
        try:
            n_slices = img.header.get_n_slices()
        except:
            pass
        else:
            n_slices = img.shape[2]
        finally:
            return n_slices

    @staticmethod
    def _get_time_repetition(in_file):
        tr = nib.load(in_file).header.get_dim_info()[2]
        if tr < 1.0 or tr > 5.0:
            raise ValueError('Aborting: strange Repeat Time (TR), got {}.'.format(tr))
        if  tr > 5.0:
            raise ValueError('Long TR often used with sparse imaging: if this is a sparse design please set the TR manually.')
        if  tr < 0.5:
            raise ValueError('Short TR may be due to DICOM-to-NIfTI conversion. Perhaps use dcm2nii.')

        return tr

    @staticmethod
    def _get_time_acquisition(in_file, TR, n_slices):
        if n_slices <= 0: # not only to avoid division by zero
            raise ValueError('Null number of slices when calculating time '
                             'acquisition for {}, got {}'.format(in_file, n_slices))

        return (TR/n_slices) * (n_slices-1)

    @staticmethod
    def _get_ref_slice(in_file, slice_order):
        if not slice_order:
            raise ValueError('Expected a list of integers as `slice_order`, got {}.'.format(slice_order))

        return slice_order[0]

    @staticmethod
    def _get_slice_order(in_file, n_slices, slice_mode):

        def read_slice_mode_byte(in_file):
            try:
                with open(in_file) as f:
                    f.seek(122)
                    slice_mode = f.read(1)
            except:
                return -1
            else:
                return slice_mode

        def get_nii_slice_times(img):
            # try get the slice times
            try:
                times = img.header.get_slice_times()
            except Exception:
                pass
            else:
                return times

        def order_from_times(times):
            return np.argsort(times) + 1

        def calculate_slice_order(n_slices, slice_mode):
            """

            Parameters
            ----------
            n_slices: int

            slice_mode: int or str
                #  0: 'unknown' : ask for automatic detection of the slice order
                #  1: 'seq_inc' : sequential ascending kNIFTI_SLICE_SEQ_INC = 1; %1,2,3,4
                #  2: 'seq_dec' : sequential descending kNIFTI_SLICE_SEQ_DEC = 2; %4,3,2,1
                #  3: 'alt_inc' : Siemens: interleaved ascending with odd number of slices, interleaved for other vendors kNIFTI_SLICE_ALT_INC = 3; %1,3,2,4
                #  4: 'alt_dec' : descending interleaved kNIFTI_SLICE_ALT_DEC = 4; %4,2,3,1
                #  5: 'alt_inc2': Siemens interleaved ascending with even number of slices kNIFTI_SLICE_ALT_INC2 = 5; %2,4,1,3
                #  6: 'alt_dec2': Siemens interleaved descending with even number of slices kNIFTI_SLICE_ALT_DEC2 = 6; %3,1,4,2

            Returns
            -------
            slice_order: list of int
            """
            mode_int = { 0: 'unknown' ,
                         1: 'seq_inc' ,
                         2: 'seq_dec' ,
                         3: 'alt_inc' ,
                         4: 'alt_dec' ,
                         5: 'alt_inc2',
                         6: 'alt_dec2',}

            if isinstance(slice_mode, int):
                slice_mode = mode_int[slice_mode]

            choices = tuple(mode_int.values())
            if slice_mode not in choices:
                raise ValueError('Expected `slice_mode` to be in {}, got {}.'.format(choices,
                                                                                     slice_mode))

            is_siemens = False
            if slice_mode in ('alt_inc2', 'alt_dec2'):
                is_siemens = True

            if 'seq' in slice_mode: # sequential
                slice_order = list(range(n_slices))
            else: # interleaved
                if is_siemens and '2' in slice_mode: #siemens and even number of slices
                    slice_order = list(range(1, n_slices, 2)) + list(range(0, n_slices, 2))
                else:
                    slice_order = list(range(0, n_slices, 2)) + list(range(0, n_slices, 2))

            if 'dec' in slice_mode: # descending
                slice_order = [n_slices - 1 - i for i in slice_order]

            return slice_order

        if slice_mode == 'unknown':
            # check if the slice times are in the NifTI header
            img = nib.load(in_file)
            times = get_nii_slice_times(img)
            if times is not None:
                return order_from_times(times)

        # read the code from the file
        if slice_mode == 'unknown':
            slice_mode = read_slice_mode_byte(in_file)

        if slice_mode > 0:
            return calculate_slice_order(n_slices, slice_mode)
        else:
            raise AttributeError("Don't have enough information to calculate the "
                                 "slice order from {}.".format(in_file))

    def _check_all_equal(self, func, error_msg=None, **kwargs):
        in_files = self._check_in_files()

        values = [func(f,  **kwargs) for f in in_files]
        if not check_equal(values):
            if error_msg is None:
                error_msg = 'The values from {} are not the same, got {{}}.'.format(func.__name__)
            raise ValueError(error_msg.format(values))

        return values[0]

    def set_time_acquisition(self):
        if self.time_acquisition is not None:
            return self.time_acquisition

        n_slices = self.set_num_slices()
        tr = self.set_time_repetition()

        error_msg = 'The time acquisition calculated from all the `in_files` are not the same, got {}.'
        self.time_acquisition = self._check_all_equal(self._get_time_acquisition,
                                                      error_msg,
                                                      TR=tr,
                                                      n_slices=n_slices)
        return self.time_acquisition

    def set_num_slices(self):
        if self.num_slices is not None:
            return self.num_slices

        error_msg = 'The number of z slices for all `in_files` are not the same, got {}.'
        self.num_slices = self._check_all_equal(self._get_n_slices, error_msg)
        return self.num_slices

    def set_time_repetition(self):
        if self.time_repetition is not None:
            return self.time_repetition

        error_msg = 'The TR calculated from all the `in_files` are not the same, got {}.'
        self.time_repetition = self._check_all_equal(self._get_time_repetition, error_msg)
        return self.time_repetition

    def set_slice_order(self):
        if self.slice_order is not None:
            return self.slice_order

        n_slices = self.set_num_slices()

        error_msg = 'The slice order for all `in_files` are not the same, got {}.'
        self.slice_order = self._check_all_equal(self._get_slice_order,
                                                 error_msg,
                                                 n_slices=n_slices,
                                                 slice_mode=self.slice_mode)

        return self.slice_order

    def set_ref_slice(self):
        """set slice order to the first slice http://www.alivelearn.net/?p=1037
           slice_order[0] is first acquired slice, set it as the refernece
        """
        if self.ref_slice is not None:
            return self.slice_order

        slice_order = self.set_slice_order()

        error_msg = 'The reference slice for all `in_files` are not the same, got {}.'
        self.ref_slice = self._check_all_equal(self._get_ref_slice, error_msg,
                                               slice_order=slice_order)
        return self.ref_slice


def slice_timing_params():
    """ Return a nipype interface to SiemensSTCParameters.

    Nipype Inputs
    -------------
    ## Mandatory:
    in_files: str or list of str
        Path to the input file(s).

    ## Optional:
    num_slices: int
        Number of slices of `in_files`.

    time_repetition: int or str
        The time repetition (TR) of the input dataset in seconds
        Default: 0
        If left to default will read the TR from the nifti image header.

    time_acquisition: int
        Time of volume acquisition. usually calculated as TR-(TR/num_slices)

    slice_order: list of int
        List of integers with the order in which slices are acquired

    ref_slice: int
        Index of the reference slice

    Nipype Outputs
    --------------
    num_slices: int
        Number of slices of `in_files`.

    ref_slice: int
        Index of the reference slice

    slice_order: list of int
        List of integers with the order in which slices are acquired

    time_acquisition: int
        Time of volume acquisition. usually calculated as TR-(TR/num_slices)

    time_repetition: int or str
        The time repetition (TR) of the input dataset in seconds
        Default: 0
        If left to default will read the TR from the nifti image header.

    Returns
    -------
    stc_iface: nipype Function

    """
    stc_pars = STCParameters()
    input_names  = [ 'in_files',
                     'num_slices',
                     'slice_order',
                     'time_repetition',
                     'time_acquisition',
                     'ref_slice',
                     'slice_mode',
                    ]

    output_names  = ['num_slices',
                     'ref_slice',
                     'slice_order',
                     'time_acquisition',
                     'time_repetition',
                    ]

    imports = ['import os.path as op',
               'import nibabel as nib',
               'import numpy as np',
               'from pypes._utils import check_equal',
              ]

    return Function(input_names=input_names,
                    output_names=output_names,
                    imports=imports,
                    function=stc_pars)
