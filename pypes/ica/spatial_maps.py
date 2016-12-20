# -*- coding: utf-8 -*-
"""
Functions to help comparing and dealing with spatial maps.
"""
import itertools

import numpy as np
import nilearn.masking as nimask
import nilearn.image as niimg
from   sklearn.metrics.pairwise import pairwise_distances
from   boyle.nifti.utils import nifti_out, thr_img, icc_img_to_zscore


@nifti_out
def spatial_map(icc, thr, mode='+'):
    """ Return the thresholded z-scored `icc`. """
    return thr_img(icc_img_to_zscore(icc), thr=thr, mode=mode).get_data()


def spatial_maps_pairwise_similarity(imgs1, imgs2, mask_file, distance='correlation'):
    """ Similarity values of each image in `imgs1` to each image in `imgs2`, both masked by `mask_file`.
    These values are based on distance metrics, specified by `distance` argument.
    The resulting similarity value is the complementary value of the distance,
    i.e., '1 - <distance value>'.
    The images in `imgs1` will be resampled to `imgs2` if their affine matrix don't match.

    Parameters
    ----------
    imgs1: list of niimg-like or 4D niimg-like

    imgs2: list of niimg-like or 4D niimg-like

    mask_file: niimg-like

    distance: str
        Valid values for `distance` are:
        From scikit-learn: ['cityblock', 'cosine', 'euclidean', 'l1', 'l2', 'manhattan'].
        From scipy.spatial.distance: ['braycurtis', 'canberra', 'chebyshev', 'correlation', 'dice', 'hamming',
                                      'jaccard', 'kulsinski', 'mahalanobis', 'matching', 'minkowski',
                                      'rogerstanimoto', 'russellrao', 'seuclidean', 'sokalmichener', 'sokalsneath',
                                      'sqeuclidean', 'yule']
                                      See the documentation for scipy.spatial.distance for details on these metrics.

    Returns
    -------
    corrs: np.ndarray
        A matrix of shape MxN, where M is len(imgs1) and N is len(imgs2).
        It contains the similarity values.
    """
    img1_ = niimg.load_img(imgs1)
    img2_ = niimg.load_img(imgs2)

    n_imgs1 = img1_.shape[-1]
    n_imgs2 = img2_.shape[-1]
    corrs = np.zeros((n_imgs1, n_imgs2), dtype=float)

    mask_trnsf = niimg.resample_to_img(mask_file, niimg.index_img(img2_, 0),
                                       interpolation='nearest',
                                       copy=True)

    for idx1, img1 in enumerate(niimg.iter_img(img1_)):
        img1_resamp = niimg.resample_to_img(img1, niimg.index_img(img2_, 0), copy=True)
        img1_masked = nimask.apply_mask(img1_resamp, mask_trnsf)

        for idx2, img2 in enumerate(niimg.iter_img(img2_)):
            img2_masked = nimask.apply_mask(img2, mask_trnsf)
            dist = pairwise_distances(img1_masked.reshape(1, -1),
                                      img2_masked.reshape(1, -1),
                                      metric=distance)

            # since this is a scalar value
            dist = dist[0][0]

            # since this is a distance, not a similarity value
            corr = 1 - dist

            # store it
            corrs[idx1, idx2] = corr

    return corrs


def spatial_maps_goodness_of_fit(rsn_imgs, spatial_maps, mask_file, rsn_thr=4.0):
    """ Goodness-of-fit values described as in Zhou et al., 2010, Brain.

    Parameters
    ----------
    rsn_imgs: list of niimg-like or 4D niimg-like
        The RSN maps. They should be thresholded beforehand if `rsn_thr` is lower or equal than 0.

    spatial_maps: list of niimg-like or 4D niimg-like

    mask_file: niimg-like
        An extra mask to apply to the thresholded RSN masks.
        This is used to exclude values outside of the RSN blobs.
        It is recommended to use a brain mask for this.

    rsn_thr: float, optional
        The threshold to apply to `rsn_imgs` to create the RSN masks.
        If rsn_thr <= 0, no thresholding will be applied.

    Returns
    -------
    gof_df: np.ndarray
        A matrix of shape MxN, where M is len(rsn_imgs) and N is len(spatial_maps).
        It contains the goodness-of-fit values.

    Notes
    -----
    "These ICN templates were thresholded at a z-score 4.0 to be visually comparable to the
    consistent ICNs published by Damoiseaux et al. (2006). A minor modification of previous
    goodness-of-fit methods (Seeley et al., 2007b, 2009) was included here for template
    matching, with goodness-of-fit scores calculated by multiplying
    (i) the average z-score difference between voxels falling within the template and
    voxels falling outside the template; and
    (ii) the difference in the percentage of positive z-score voxels inside and outside the template.

    This goodness-of-fit algorithm proves less vulnerable to inter-subject variability in shape,
    size, location and strength of each ICN", than the one published in Seeley et al., 2007b.
    This latter method only uses the value in (i).

    Extracted from Zhou et al., 2010, Brain.
    """
    rsn_img = niimg.load_img(rsn_imgs)
    spm_img = niimg.load_img(spatial_maps)

    n_rsns = rsn_img.shape[-1]
    n_ics  =  spm_img.shape[-1]
    gofs   = np.zeros((n_rsns, n_ics), dtype=float)

    # threshold the RSN templates
    if rsn_thr > 0:
        thr_rsns = (spatial_map(rsn, thr=rsn_thr, mode='+-')
                    for rsn in niimg.iter_img(rsn_img))
    else:
        thr_rsns = rsn_img

    # for each RSN template and IC image
    iter_rsn_ic = itertools.product(enumerate(niimg.iter_img(thr_rsns)),
                                    enumerate(niimg.iter_img( spm_img)))

    for (rsn_idx, rsn), (ic_idx, ic) in iter_rsn_ic:

        ref_vol = rsn.get_data()
        rsn_vol = np.zeros(rsn.shape, dtype=int)

        #rsn_in  = niimg.math_img('np.abs(img) > 0', img=rsn)
        rsn_in = rsn_vol.copy()
        rsn_in[np.abs(ref_vol) > 0] = 1

        #rsn_out = niimg.math_img('img == 0', img=rsn)
        rsn_out = rsn_vol.copy()
        rsn_out[ref_vol == 0] = 1

        if mask_file is not None:
            # rsn_out = niimg.math_img('mask * img', mask=rsn_brain_mask, img=rsn_out)
            rsn_brain_mask = niimg.resample_to_img(mask_file, rsn,
                                                   interpolation='nearest')
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
        gof_term1 = zscore_in.mean() - zscore_out.mean()

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


from sklearn.linear_model import Ridge


def _compute_loadings(components, data):
    """
    DONT USE THIS YET!

    Parameters
    ----------
    components:

    data:

    Returns
    -------

    """
    ridge = Ridge(fit_intercept=None, alpha=1e-8)
    ridge.fit(components.T, np.asarray(data.T))
    loadings = ridge.coef_.T

    # unit-variance scaling
    S = np.sqrt(np.sum(loadings ** 2, axis=0))
    S[S == 0] = 1
    loadings /= S[np.newaxis, :]
    return loadings


def spatiotemporal_regression(components, imgs, mask, confounds=None):
    """Project the data into a reduced representation

    DONT USE THIS YET!

    Parameters
    ----------
    components:

    imgs: iterable of Niimg-like objects
        See http://nilearn.github.io/manipulating_images/input_output.html.
        Data to be projected

    confounds: CSV file path or 2D matrix
        This parameter is passed to nilearn.signal.clean. Please see the
        related documentation for details

    Returns
    ----------
    loadings: list of 2D ndarray,
        For each subject, each sample, loadings for each decomposition
        components
        shape: number of subjects * (number of scans, number of regions)

    """
    components_img_ = self.masker_.inverse_transform(self.components_)

    nifti_maps_masker = NiftiMapsMasker(
        components_img_, self.masker_.mask_img_,
        resampling_target='maps')

    nifti_maps_masker.fit()

    if confounds is None:
        confounds = [None] * len(imgs)

    return [nifti_maps_masker.transform(img, confounds=confound)
            for img, confound in zip(imgs, confounds)]


def inverse_transform(self, loadings):
    """Use provided loadings to compute corresponding linear component
    combination in whole-brain voxel space

    DONT USE THIS YET

    Parameters
    ----------
    loadings: list of numpy array (n_samples x n_components)
        Component signals to tranform back into voxel signals

    Returns
    -------
    reconstructed_imgs: list of nibabel.Nifti1Image
        For each loading, reconstructed Nifti1Image

    """
    if not hasattr(self, 'components_'):
        ValueError('Object has no components_ attribute. This is either '
                   'because fit has not been called or because'
                   '_DecompositionEstimator has directly been used')
    self._check_components_()
    components_img_ = self.masker_.inverse_transform(self.components_)
    nifti_maps_masker = NiftiMapsMasker(
        components_img_, self.masker_.mask_img_,
        resampling_target='maps')
    nifti_maps_masker.fit()
    # XXX: dealing properly with 2D/ list of 2D data?
    return [nifti_maps_masker.inverse_transform(loading)
            for loading in loadings]

