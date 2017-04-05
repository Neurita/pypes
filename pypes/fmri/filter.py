# -*- coding: utf-8 -*-
"""
fMRI timeseries filtering helpers.
"""


def bandpass_filter(files, lowpass_freq=0.1, highpass_freq=0.01, tr=2):
    """Bandpass filter the input files

    Parameters
    ----------
    files: list of str
        List 4d nifti file paths.

    lowpass_freq: float
        Cutoff frequency for the low pass filter (in Hz).

    highpass_freq: float
        Cutoff frequency for the high pass filter (in Hz).

    tr: float
        The repetition time in seconds. The inverse of sampling rate (in Hz).
    """
    import os

    import nibabel as nb
    import numpy as np
    from   nipype.utils.filemanip import (filename_to_list,
                                          list_to_filename,
                                          split_filename)

    fs = 1./tr

    out_files = []
    for filename in filename_to_list(files):
        path, name, ext = split_filename(filename)
        out_file = os.path.join(os.getcwd(), name + '_bandpassed' + ext)

        img = nb.load(filename)
        timepoints = img.shape[-1]
        F = np.zeros((timepoints))

        lowidx = int(timepoints / 2) + 1
        if lowpass_freq > 0:
            lowidx = np.round(float(lowpass_freq) / fs * timepoints)

        highidx = 0
        if highpass_freq > 0:
            highidx = np.round(float(highpass_freq) / fs * timepoints)
        F[int(highidx):int(lowidx)] = 1
        F = ((F + F[::-1]) > 0).astype(int)
        data = img.get_data()
        if np.all(F == 1):
            filtered_data = data
        else:
            filtered_data = np.real(np.fft.ifftn(np.fft.fftn(data) * F))
        img_out = nb.Nifti1Image(filtered_data, img.affine, img.header)
        img_out.to_filename(out_file)
        out_files.append(out_file)

    return list_to_filename(out_files)
