# -*- coding: utf-8 -*-
"""
Helper functions to read, threshold, and build IC loadings results.
"""
import pandas as pd
import nilearn.image as niimg
from   nilearn.image import iter_img
from   boyle.nifti.roi import largest_connected_component
from   boyle.nifti.utils import filter_icc


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


def build_raw_loadings_table(loads, patids):
    """ Build a spreadsheet-ready pandas.DataFrame with the content of the
    loadings matrix and subjects ids.
    """
    loadings = []
    for p in range(len(patids)):
        patid = patids[p]
        loadings.append([patid] + list(loads[p, :]))

    # set the column names
    n_cols = loads.shape[1]
    cols = ['subject_id'] + list(range(1, n_cols+1))

    # fill the df
    return pd.DataFrame.from_records(loadings, columns=cols)


def add_groups_to_loadings_table(df, groups):
    """ Merge `df` and `groups` on 'subject_id' and sort by 'group'.
    Note: `groups` must have 'subject_id' (as well as `df`) and 'group' columns.
    """
    if groups is not None:
        df = pd.merge(df, groups, on='subject_id')
        df = df.sort_values(by=['group', 'subject_id'], ascending=[True, True])

    df.loc[:, 'group'] = df.loc[:, 'group'].astype('category')
    return df


def filter_ics(comps_img, mask, zscore=2., mode='+-'):
    """
    Generator for masking and thresholding each IC spatial map.

    Parameters
    ----------
    comps_img: img-like
        The 'raw' ICC maps image.

    mask: img-like
        If not None. Will apply this masks in the end of the process.

    thr: float
        The threshold value.

    zscore: bool
        If True will calculate the z-score of the ICC before thresholding.

    mode: str
        Choices: '+' for positive threshold,
                 '+-' for positive and negative threshold and
                 '-' for negative threshold.

    Returns
    -------
    icc_filts: list of nibabel.NiftiImage
        Thresholded and masked ICCs.
    """
    # store the average value of the blob in a list
    mask = niimg.load_img(mask)
    for i, icimg in enumerate(iter_img(comps_img)):
        yield filter_icc(icimg, mask=mask, thr=zscore, zscore=True, mode=mode)
