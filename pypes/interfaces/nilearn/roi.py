# -*- coding: utf-8 -*-
"""
Utilities to use as mask and extract labels and ROIs from images.
"""
import nilearn.image as niimg
import numpy as np


def spread_labels(labels_img, roi_values=None, background_label=0):
    """ Spread each ROI in labels_img into a 4D image.

    Parameters
    ----------
    labels_img: 3D Niimg-like object
        See http://nilearn.github.io/manipulating_images/input_output.html. Region definitions, as one image of labels.

    roi_values: List of int or the values in the labels_img, optional
        The values of `labels_img` that you want to use.
        If None, will extract all the values except from `background_label`.

    background_label: number, optional
        Label used in labels_img to represent background.

    Return
    ------
    4d_labels_img: 4D Niimg-like object
    """
    lbls_img = niimg.load_img(labels_img)
    lbls_vol = lbls_img.get_data()

    if roi_values is None:
        roi_values = np.unique(lbls_vol)
        roi_values = roi_values[roi_values != background_label]

    atlas_roi_vols = [(lbls_vol == val) * val for val in roi_values]

    atlas_roi_imgs = [niimg.new_img_like(lbls_img, vol) for vol in atlas_roi_vols]
    # atlas_roi_imgs = [niimg.new_img_like(lbls_img, (lbls_vol == val) * val) for val in roi_values]

    return niimg.concat_imgs(atlas_roi_imgs)