# -*- coding: utf-8 -*-
"""
Nipype workflows to process anatomical MRI.
"""
import json
from collections import OrderedDict
from typing import Iterable

import pandas as pd
from nipype.interfaces.nipy.preprocess import Trim
from nipype.utils.filemanip import filename_to_list
from pydicom import read_file

import hansel


def get_info(dicom_files):
    """Given a Siemens dicom file return metadata

    Returns
    -------
    RepetitionTime
    Slice Acquisition Times
    Spacing between slices
    """
    try:
        from dcmstack.extract import default_extractor
    except:
        raise ImportError("To use this function install dcmstack: pip install "
                          "git+https://github.com/moloney/dcmstack.git@c12d27d2c802d75a33ad70110124500a83e851ee#egg=dcmstack")

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


def motion_stats_sheet(motion_filepath: [str, hansel.Crumb], crumb_fields: Iterable[str]):
    """ Return a pandas.DataFrame with some of the motion statistics obtained from the
    `statistics_files` output of the nipype.RapidArt found in the hansel.Crumb `motion_file_cr`.

    Parameters
    ----------
    motion_filepath
    crumb_fields:

    Returns
    -------
    df: pandas.DataFrame

    Examples
    --------
    >>> motion_stats_sheet(motion_file_cr="~/data/out/{group}/{sid}/session_0/rest/artifact_stats/motion_stats.json", \
    >>>                    crumb_fields=['group', 'patient_id', 'session'])
    """
    def get_motion_record(mtn_file_cr, crumb_fields):
        """ Return an OrderedDict of the information found in the `mtn_file_cr` and also
        `crumb_fields` Crumb argument values."""
        stats = json.load(open(str(mtn_file_cr)))

        outliers = stats[1]
        motion_norm = stats[3]['motion_norm']

        #outliers_hdr = list(outliers.keys())
        motion_hdr   = ['{}_motion_norm'.format(k) for k in motion_norm.keys()]

        mtn_record = OrderedDict()
        for fn in crumb_fields:
            mtn_record[fn] = mtn_file_cr[fn][0]
        mtn_record.update(outliers)

        for hdr, fn in zip(motion_hdr, motion_norm):
            mtn_record[hdr] = motion_norm[fn]

        return mtn_record

    # process the input
    if isinstance(motion_filepath, str):
        motion_filepath = hansel.Crumb(motion_filepath)

    # create the motion records
    motionstats = [get_motion_record(stats_file, crumb_fields) for stats_file in motion_filepath.ls()]

    # create a pandas Dataframe out of it
    df = pd.DataFrame.from_records(motionstats, columns=motionstats[0].keys())

    # return the dataframe
    return df
