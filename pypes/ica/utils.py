# -*- coding: utf-8 -*-
"""
Helper functions to read, threshold, and build IC loadings results.
"""
import nilearn.image   as niimg
from   nilearn.image   import iter_img
from   boyle.nifti.roi import largest_connected_component


def get_largest_blobs(ic_maps):
    """ Generator for the largest blobs in each IC spatial map.
    These should be masked and thresholded.

    Parameters
    ----------
    ic_maps: sequence of niimg-like

    Returns
    -------
    blobs: generator of niimg-like
    """
    # store the average value of the blob in a list
    for i, icimg in enumerate(iter_img(ic_maps)):
        yield niimg.new_img_like(icimg, largest_connected_component(icimg.get_data()))
