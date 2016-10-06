# -*- coding: utf-8 -*-
"""
Resting state networks post-hoc analysis functions.
"""
from collections import OrderedDict
import itertools

import numpy as np
import pandas as pd
import nilearn.plotting as niplot
import nilearn.masking as nimask
import nilearn.image as niimg
from   sklearn.metrics.pairwise import pairwise_distances
from   boyle.nifti.utils import spatial_map


class RestingStateNetworks:
    """ A container class to parse and return useful values/images
    from RSN templates.

    Parameters
    ----------
    img_file: str
        Path to a img-like RSN templates file.

    txt_file: str
        Path to a text file with the labels for `img_file`.
        The expected structure of each line of 
        this file is  <name of the RSN>, <list of indices>.
        For example:
        
        Basal Ganglia, 21
        Auditory, 17
        Sensorimotor, 7, 23, 24, 38, 56, 29
        Visual, 46, 64, 67, 48, 39, 59
        Default-Mode, 50, 53, 25, 68
        Attentional, 34, 60, 52, 72, 71, 55
        Frontal, 42, 20, 47, 49

    start_from_one: bool
        If True it means that the `txt_file` volume indices start from 1.
        From 0 otherwise. Be careful, the default is True!
    """
    def __init__(self, img_file, txt_file, start_from_one=True):
        self._img_file       = img_file
        self._txt_file       = txt_file
        self._start_from_one = start_from_one
        
        self.network_names = self._network_names()
        self._img          = niimg.load_img(self._img_file)
        self._self_check()

    def iter_networks(self):
        """Yield idx (what is in the text_file) and
        image of each RSN volume from the text file."""
        for idx, name in self.network_names.items():
            yield idx, self._get_img(idx)

    def _network_names(self):
        """Return OrderedDict[int]->str, with the index and the name of each
        RSN."""
        names_idx = self._read_labels_file()
        return OrderedDict([(idx, name) for name, idxs in names_idx.items()
                                        for idx in idxs])

    def _get_img(self, network_index):
        """ Return one RSN given the index in the labels file."""
        img_idx = self._img_index(network_index)
        return niimg.index_img(self._img, img_idx)

    def _img_index(self, network_index):
        """Return the correspoding image index for the given network index."""
        if self._start_from_one:
            return network_index - 1

        return network_index

    def _self_check(self):
        """Simple content check."""
        n_labels = len(self.network_names)
        n_images = self._img.shape[-1]

        if n_labels == n_images:
            return

        if n_labels > n_images:
            raise ValueError('The number of labels is larger than the number '
                             'of images. Got {} and {}.'.format(n_labels, n_images))

        # print('The number of volumes in the image is different from the number '
        #       'of labels in the text file.\n I am going to use only the ones '
        #       ' in the text file.')

    def _read_labels_file(self):
        """ Read the text file and return a dict[str->List[int]] with network
        names and blob indices.
        """
        lines = [l.rstrip('\n') for l in open(self._txt_file).readlines()]
        
        netblobs = OrderedDict()
        for l in lines:
            pcs = l.split(',')
            netname = pcs[0]
            blobs   = [int(idx) for idx in pcs[1:]]

            netblobs[netname] = blobs
        return netblobs

    def plot_all(self):
        names = self.network_names
        for idx, rsn in enumerate(niimg.iter_img(self._img)):
            disp = niplot.plot_roi(rsn, title=names.get(idx, None))

    def join_networks(self, network_indices):
        """Return a NiftiImage containing a binarised version of the sum of
        the RSN images of each of the `network_indices`."""
        oimg = self._get_img(network_indices[0]).get_data()
        for idx in network_indices[1:]:
            oimg += self._get_img(idx).get_data()

        return niimg.new_img_like(self._get_img(network_indices[0]),
                                  oimg.astype(bool))

    def __iter__(self):
        return (img for idx, img in self.iter_networks())

    def __len__(self):
        return len(self.network_names)


def spatial_maps_pairwise_similarity(rsn_imgs, ic_imgs, mask_file, distance='correlation'):
    """ Similarity values of each RSN to each IC map in `ic_imgs` masked by `mask_file`.
    These values are based on distance metrics, specified by `distance` argument.
    The resulting similarity value is the complementary value of the distance,
    i.e., '1 - <distance value>'.

    Parameters
    ----------
    rsn_imgs: list of niimg-like or 4D niimg-like

    ic_imgs: list of niimg-like or 4D niimg-like

    mask_file: niimg-like

    distance: str
        Valid values for metric are:
        From scikit-learn: ['cityblock', 'cosine', 'euclidean', 'l1', 'l2', 'manhattan'].
        From scipy.spatial.distance: ['braycurtis', 'canberra', 'chebyshev', 'correlation', 'dice', 'hamming',
                                      'jaccard', 'kulsinski', 'mahalanobis', 'matching', 'minkowski',
                                      'rogerstanimoto', 'russellrao', 'seuclidean', 'sokalmichener', 'sokalsneath',
                                      'sqeuclidean', 'yule']
                                      See the documentation for scipy.spatial.distance for details on these metrics.

    Returns
    -------
    corrs: np.ndarray
        A matrix of shape MxN, where M is len(rsn_imgs) and N is len(ic_imgs).
        It contains the correlation values.
    """
    rsn_img = niimg.load_img(rsn_imgs)
    ic_img  = niimg.load_img( ic_imgs)

    n_rsns = rsn_img.shape[-1]
    n_ics  =  ic_img.shape[-1]
    corrs = np.zeros((n_rsns, n_ics), dtype=float)

    mask_trnsf = niimg.resample_to_img(mask_file, niimg.index_img(ic_img, 0),
                                       interpolation='nearest',
                                       copy=True)

    for rsn_idx, rsn in enumerate(niimg.iter_img(rsn_img)):
        rsn_transf = niimg.resample_to_img(rsn, niimg.index_img(ic_img, 0), copy=True)
        rsn_masked = nimask.apply_mask(rsn_transf, mask_trnsf)

        for ic_idx, ic in enumerate(niimg.iter_img(ic_img)):
            ic_masked = nimask.apply_mask(ic, mask_trnsf)
            dist = pairwise_distances(rsn_masked.reshape(1, -1),
                                       ic_masked.reshape(1, -1),
                                      metric=distance)

            # since this is a scalar value
            dist = dist[0][0]

            # since this is a distance based on correlation, not a correlation value
            corr = 1 - dist

            # store it
            corrs[rsn_idx, ic_idx] = corr

    return corrs


def spatial_maps_goodness_of_fit(rsn_imgs, ic_imgs, mask_file, rsn_thr=4.0):
    """ Goodness-of-fit values described as in Zhou et al., 2010, Brain.

    Parameters
    ----------
    rsn_imgs: list of niimg-like or 4D niimg-like

    ic_imgs: list of niimg-like or 4D niimg-like

    mask_file: niimg-like

    rsn_thr: float

    Returns
    -------
    gof_df: np.ndarray
        A matrix of shape MxN, where M is len(rsn_imgs) and N is len(ic_imgs).
        It contains the goodness-of-fit values.

    Notes
    -----
    These ICN templates were thresholded at a z-score 4.0 to be visually comparable to the
    consistent ICNs published by Damoiseaux et al. (2006). A minor modification of previous
    goodness-of-fit methods (Seeley et al., 2007b, 2009) was included here for template
    matching, with goodness-of-fit scores calculated by multiplying
    (i) the average z-score difference between voxels falling within the template and
    voxels falling outside the template; and
    (ii) the difference in the percentage of positive z-score voxels inside and outside the template.

    This goodness-of-fit algorithm proves less vulnerable to inter-subject variability in shape,
    size, location and strength of each ICN.
    """
    rsn_img = niimg.load_img(rsn_imgs)
    ic_img  = niimg.load_img( ic_imgs)

    n_rsns = rsn_img.shape[-1]
    n_ics  =  ic_img.shape[-1]
    gofs   = np.zeros((n_rsns, n_ics), dtype=float)

    # threshold the RSN templates
    if rsn_thr > 0:
        thr_rsns = (spatial_map(rsn, thr=rsn_thr, mode='+-')
                    for rsn in niimg.iter_img(rsn_img))
    else:
        thr_rsns = rsn_img

    # for each RSN template and IC image
    iter_rsn_ic = itertools.product(enumerate(niimg.iter_img(thr_rsns)),
                                    enumerate(niimg.iter_img( ic_img )))

    for (rsn_idx, rsn), (ic_idx, ic) in iter_rsn_ic:
        # prepare the RSN masks
        rsn_brain_mask = niimg.resample_to_img(mask_file, rsn,
                                               interpolation='nearest')

        ref_vol = rsn.get_data()
        rsn_vol = np.zeros(rsn.shape, dtype=int)
        #rsn_in  = niimg.math_img('np.abs(img) > 0', img=rsn)
        rsn_in = rsn_vol.copy()
        rsn_in[np.abs(ref_vol) > 0] = 1

        #rsn_out = niimg.math_img('img == 0', img=rsn)
        rsn_out = rsn_vol.copy()
        rsn_out[ref_vol == 0] = 1

        #rsn_out = niimg.math_img('mask * img', mask=rsn_brain_mask, img=rsn_out)
        rsn_out = rsn_brain_mask.get_data() * rsn_out

        # convert the mask arrays to image in order to resample
        rsn_in  = niimg.new_img_like(rsn, rsn_in)
        rsn_out = niimg.new_img_like(rsn, rsn_out)

        rsn_in  = niimg.resample_to_img(rsn_in,  ic, interpolation='nearest')
        rsn_out = niimg.resample_to_img(rsn_out, ic, interpolation='nearest')

        # apply the mask
        #zscore_in  = niimg.math_img('mask * img', mask=rsn_in,  img=ic).get_data()
        zscore_in = rsn_in.get_data() * ic.get_data()

        #zscore_out = niimg.math_img('mask * img', mask=rsn_out, img=ic).get_data()
        zscore_out = rsn_out.get_data() * ic.get_data()

        #gof_term1
        # calculate the the average z-score difference between voxels falling
        # within the template and voxels falling outside the template
        gof_term1  = zscore_in.mean() - zscore_out.mean()

        #gof_term2
        # the difference in the percentage of positive z-score voxels inside and outside the template.
        n_pos_zscore_in  = np.sum(zscore_in  > 0)
        n_pos_zscore_out = np.sum(zscore_out > 0)
        n_pos_zscore_tot = n_pos_zscore_in + n_pos_zscore_out

        if n_pos_zscore_tot != 0:
            n_pos_zscore_pcnt = 100 / n_pos_zscore_tot
            gof_term2 = (n_pos_zscore_in - n_pos_zscore_out) * n_pos_zscore_pcnt
        else:
            gof_term2 = 0

        # global gof
        gof = gof_term1 * gof_term2

        # add the result
        gofs[rsn_idx][ic_idx] = gof

    return gofs


def label_pairwise_measure(values, axis0_labels, axis1_labels, **kwargs):
    """ Return a pandas.DataFrame with the content of `values`.
    This DataFrame will have each row labeled by `axis0_labels`, and
    each column by `axis1_labels`.

    Parameters
    ----------
    values: np.ndarray
        MxN matrix
        M is the length of axis0_labels.
        N is the lenght of axis1_labels.

    axis0_labels: array or any type
        This will be used as index in the DataFrame construction.

    axis1_labels: array or any type
        This will be used as column in the DataFrame construction.

    kwargs: keyword arguments
        Additional columns to be added to the resulting DataFrame.
        The argument names are the column name and the values must be sequences of length M.

    Returns
    -------
    labeled_measure: pandas.DataFrame

    Notes
    -----
    My previous solution was more complex and I had this function
    prepared before realizing this was ridiculously simple.
    """
    df = pd.DataFrame(values, index=axis0_labels, columns=axis1_labels)

    for k, v in kwargs.items():
        if len(v) != len(axis0_labels):
            raise AttributeError('The value for argument {} should have length {} but has '
                                 'length {}.'.format(k, len(axis0_labels), len(v)))
        df[k] = v

    return df


def nd_vector_correlations(data, vector, n=4):
    """ Calculate the pairwise correlation between `vector` and each element
    of the last dimension of `data`.

    Parameters
    ----------
    data: numpy.ndarray
        nD array

    vector: numpy.array
        1D vector

    n: int

    Returns
    -------
    corrs: numpy.array
        (n-1)D array
    """
    assert data.ndim == n
    assert data.shape[-1] == list(vector.shape)[0]

    trid_shape = data.shape[:n-1]

    n_voxels = np.prod(trid_shape)
    n_ts     = data.shape[-1]

    fourthd_vecs = data.reshape(n_voxels, n_ts)
    corrs        = - pairwise_distances(fourthd_vecs, vector[np.newaxis], 'correlation') + 1
    corrs[np.isnan(corrs)] = 0

    return corrs.reshape(trid_shape)


# def subj_ic_spatial_map(fourd_img, ic):
#     """
#     """
#     img = niimg.load_img(fourd_img)
#
#     corrs = _4d_vector_correlations(img.get_data(), ic)
#
#     return niimg.new_img_like(img, corrs)


# corrs_df = label_pairwise_measure(correlations, allens.network_names.keys(),
#                                   allens.network_names.values())
