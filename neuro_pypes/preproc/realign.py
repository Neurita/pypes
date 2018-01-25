# -*- coding: utf-8 -*-
"""
Helper functions for motion correction of 4D images.
"""
from nipype.interfaces.base import traits
from nipype.interfaces.nipy import SpaceTimeRealigner


def nipy_motion_correction(in_file=traits.Undefined):
    """
    Simultaneous motion and slice timing correction algorithm
    If slice_times is not specified, this algorithm performs spatial motion correction
    This interface wraps nipy’s SpaceTimeRealign algorithm [Roche2011] or
    simply the SpatialRealign algorithm when timing info is not provided.

    http://www.mit.edu/~satra/nipype-nightly/interfaces/generated/nipype.interfaces.nipy.preprocess.html#spacetimerealigner

    Returns
    -------
    realign: nipype Interface

    Notes
    -----
    Roche A. A four-dimensional registration algorithm with application to joint correction of motion and slice
    timing in fMRI. IEEE Trans Med Imaging. 2011 Aug;30(8):1546-54. http://dx.doi.org/10.1109/TMI.2011.2131152
    """
    realigner = SpaceTimeRealigner()
    realigner.inputs.in_file = in_file
    return realigner
