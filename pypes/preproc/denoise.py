# -*- coding: utf-8 -*-
"""
Denoise and motion correction helper functions
"""
from boyle.nifti.utils import nifti_out


@nifti_out
def nlmeans_denoise_img(img, mask, N=4):
    """ Apply dipy nlmeans denoising to the img. Useful for diffusion images.
    Parameters
    ----------
    img: nibabel.Nifti1Image
        The diffusion image

    mask: nibabel.Nifti1Image
        A brain mask image.

    N: int
        Number of arrays of the head coil used to acquired the image.

    Returns
    -------
    den_img: nibabel.Nifti1Image
        A denoised nifti image object with the same headers
        and affine as `img`.
    """
    from dipy.denoise.nlmeans import nlmeans
    from dipy.denoise.noise_estimate import estimate_sigma

    data = img.get_data()
    msk  = mask.get_data()

    sigma = estimate_sigma(data, N=N)
    return nlmeans(data, sigma=sigma, mask=msk)


def reslice_img(img, new_zooms=None, order=3):
    """ Performs regridding of an image to set isotropic voxel sizes using dipy.
    If the file has already isotropic voxels, will return a copy of the same image.

    Parameters
    ----------
    img: nibabel.Nifti1Image
        The diffusion image.

    new_zooms : tuple, shape (3,)
        new voxel size for (i,j,k) after resampling

    order : int, from 0 to 5
       order of interpolation for resampling/reslicing, 0 nearest interpolation, 1 trilinear etc..
       if you don’t want any smoothing 0 is the option you need.

    Returns
    -------
    nu_img: nibabel.Nifti1Image
        A isotropic voxel version from `img`.
    """
    import numpy as np
    import nibabel as nib
    from dipy.align.reslice import reslice

    # read the data
    data = img.get_data()
    img_zooms = img.header.get_zooms()[:3]

    # check if already isotropic voxels
    all_equal = len(np.unique(img_zooms)) == 1
    if all_equal:
        return nib.Nifti1Image(img.get_data(), affine=img.affine, header=img.header)

    # set new_zooms parameter
    if new_zooms is None:
        minzoom = np.array(img_zooms).min()
        new_zooms = tuple(np.ones((3,)) * minzoom)

    # reslice it
    nu_data, nu_affine = reslice(data=data, affine=img.affine,
                                 zooms=img_zooms,
                                 new_zooms=new_zooms,
                                 order=order)

    # create the new image object
    header = img.header.copy()
    tmp_zooms = np.array(header.get_zooms())
    tmp_zooms[:3] = new_zooms[0]
    header.set_zooms(tuple(tmp_zooms))
    header.set_data_shape(nu_data.shape)
    header.set_xyzt_units('mm')
    nu_img = nib.Nifti1Image(nu_data.astype(header.get_data_dtype()), nu_affine, header)

    return nu_img


def motion_regressors(motion_params, order=0, derivatives=1):
    """Compute motion regressors upto given order and derivative

    motion + d(motion)/dt + d2(motion)/dt2 (linear + quadratic)
    """
    import os
    import numpy as np
    from nipype.utils.filemanip import filename_to_list

    out_files = []
    for idx, filename in enumerate(filename_to_list(motion_params)):
        params = np.genfromtxt(filename)
        out_params = params
        for d in range(1, derivatives + 1):
            cparams = np.vstack((np.repeat(params[0, :][None, :], d, axis=0),
                                 params))
            out_params = np.hstack((out_params, np.diff(cparams, d, axis=0)))
        out_params2 = out_params
        for i in range(2, order + 1):
            out_params2 = np.hstack((out_params2, np.power(out_params, i)))
        filename = os.path.join(os.getcwd(), "motion_regressor%02d.txt" % idx)
        np.savetxt(filename, out_params2, fmt=b"%.10f")
        out_files.append(filename)
    return out_files


def create_regressors(motion_params, comp_norm, outliers, detrend_poly=None):
    """Builds a regressor set comprising motion parameters, composite norm and
    outliers.
    The outliers are added as a single time point column for each outlier

    Parameters
    ----------
    motion_params: a text file containing motion parameters and its derivatives
    comp_norm: a text file containing the composite norm
    outliers: a text file containing 0-based outlier indices
    detrend_poly: number of polynomials to add to detrend

    Returns
    -------
    components_file: a text file containing all the regressors
    """
    import os
    import numpy as np
    from nipype.utils.filemanip import filename_to_list
    from scipy.special import legendre

    out_files = []
    for idx, filename in enumerate(filename_to_list(motion_params)):
        params = np.genfromtxt(filename)
        norm_val = np.genfromtxt(filename_to_list(comp_norm)[idx])
        out_params = np.hstack((params, norm_val[:, None]))
        try:
            outlier_val = np.genfromtxt(filename_to_list(outliers)[idx])
        except IOError:
            outlier_val = np.empty((0))
        for index in np.atleast_1d(outlier_val):
            outlier_vector = np.zeros((out_params.shape[0], 1))
            outlier_vector[index] = 1
            out_params = np.hstack((out_params, outlier_vector))
        if detrend_poly:
            timepoints = out_params.shape[0]
            X = np.empty((timepoints, 0))
            for i in range(detrend_poly):
                X = np.hstack((X, legendre(
                    i + 1)(np.linspace(-1, 1, timepoints))[:, None]))
            out_params = np.hstack((out_params, X))
        filename = os.path.join(os.getcwd(), "filter_regressor%02d.txt" % idx)
        np.savetxt(filename, out_params, fmt=b"%.10f")
        out_files.append(filename)
    return out_files


def extract_noise_components(realigned_file, mask_file, num_components=5,
                             extra_regressors=None):
    """Derive components most reflective of physiological noise.
    Parameters
    ----------
    realigned_file: a 4D Nifti file containing realigned volumes
    mask_file: a 3D Nifti file containing white matter + ventricular masks
    num_components: number of components to use for noise decomposition
    extra_regressors: additional regressors to add
    Returns
    -------
    components_file: a text file containing the noise components
    """
    import os
    import nibabel as nb
    import numpy as np
    import scipy as sp
    from   nipype.utils.filemanip import filename_to_list

    imgseries = nb.load(realigned_file)
    components = None
    for filename in filename_to_list(mask_file):
        mask = nb.load(filename).get_data()
        if len(np.nonzero(mask > 0)[0]) == 0:
            continue
        voxel_timecourses = imgseries.get_data()[mask > 0]
        voxel_timecourses[np.isnan(np.sum(voxel_timecourses, axis=1)), :] = 0
        # remove mean and normalize by variance
        # voxel_timecourses.shape == [nvoxels, time]
        X = voxel_timecourses.T
        stdX = np.std(X, axis=0)
        stdX[stdX == 0] = 1.
        stdX[np.isnan(stdX)] = 1.
        stdX[np.isinf(stdX)] = 1.
        X = (X - np.mean(X, axis=0)) / stdX
        u, _, _ = sp.linalg.svd(X, full_matrices=False)
        if components is None:
            components = u[:, :num_components]
        else:
            components = np.hstack((components, u[:, :num_components]))
    if extra_regressors:
        regressors = np.genfromtxt(extra_regressors)
        components = np.hstack((components, regressors))
    components_file = os.path.join(os.getcwd(), 'noise_components.txt')

    np.savetxt(components_file, components, fmt=b"%.10f")
    return components_file
