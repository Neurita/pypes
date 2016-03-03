# -*- coding: utf-8 -*-
"""
Nipype workflows to process anatomical MRI.
"""

from nipype.utils.filemanip import filename_to_list
from nipype.interfaces.nipy.preprocess import Trim
from nipype.interfaces.c3 import C3dAffineTool

from dicom import read_file


def get_info(dicom_files):
    from dcmstack.extract import default_extractor
    """Given a Siemens dicom file return metadata

    Returns
    -------
    RepetitionTime
    Slice Acquisition Times
    Spacing between slices
    """
    meta = default_extractor(read_file(filename_to_list(dicom_files)[0],
                                       stop_before_pixels=True,
                                       force=True))
    return (meta['RepetitionTime'] / 1000., meta['CsaImage.MosaicRefAcqTimes'],
            meta['SpacingBetweenSlices'])


def trim(begin_idx=6):
    """Simple interface to trim a few volumes from a 4d fmri nifti file.

    :param begin_idx: int
        Remove first `begin_idx` volumes.

    Notes
    -----
    http://www.mit.edu/~satra/nipype-nightly/interfaces/generated/nipype.interfaces.nipy.preprocess.html#trim
    """
    trim = Trim()
    #trim.inputs.in_file = 'functional.nii'
    trim.inputs.begin_index = begin_idx
    return trim