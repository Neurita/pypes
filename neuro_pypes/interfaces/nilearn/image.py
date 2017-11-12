# -*- coding: utf-8 -*-
"""
Utilities to use nilearn.image from nipype
"""


def ni2file(**presuffixes):
    from functools import wraps

    """ Pick the nibabel image container output from `f` and stores it in a file.
    If the shape attribute of the container is True, will save it into a file, otherwise
    will directly return the scalar value.

    To know the path where it has to save the file it looks into argument
    values in this order:
    - `out_file` kwarg in the call to `f`,
    - `out_file` kwarg in the `ni2file` decorator,
    - the first argument in the function call.
    In the last case a presuffix must be defined in the decorator to avoid overwriting
    an existing file.
    """
    def _pick_an_input_file(*args, **kwargs):
        """Assume that either the first arg or the first kwarg is an input file."""
        if args:
            return args[0]
        else:
            return list(kwargs.values())[0]

    def nifti_out(f):
        @wraps(f)
        def wrapped(*args, **kwargs):
            res_img = f(*args, **kwargs)
            if isinstance(res_img, list):
                if len(res_img) == 1:
                    res_img = res_img[0]
                else:
                    return res_img

            if not res_img.shape: # the result is a scalar value
                return res_img.get_data().flatten()[0]

            import os.path as op
            from nipype.utils.filemanip import fname_presuffix

            out_file = kwargs.get('out_file', presuffixes.pop('out_file', None))
            if out_file is not None:
                if not presuffixes and op.exists(out_file):
                    raise IOError('The file {} already exists, please add a presuffix to the'
                                  'decorator.'.format(out_file))

                out_file = fname_presuffix(out_file, **presuffixes)
            else:
                in_file = kwargs.get('in_file', None)
                if in_file is None:
                    in_file = _pick_an_input_file(*args, **kwargs)

                if not op.exists(in_file):
                    raise IOError('Expected an existing file to use as reference for'
                                  ' the output file name, got {}.'.format(in_file))

                out_file = fname_presuffix(op.basename(in_file), **presuffixes)

            if not out_file:
                raise ValueError("Could not find a output file name for this function: "
                                " {}({}, {}).".format(f.__name__, *args, **kwargs))

            res_img.to_filename(out_file)

            return op.abspath(out_file)

        return wrapped
    return nifti_out


@ni2file(out_file='nilearn_maths.nii.gz')
def math_img(formula, out_file='', **imgs):
    """ Use nilearn.image.math_img.
    This function in addition allows imgs to contain numerical scalar values.

    Returns
    -------
    out_file: str
        The absolute path to the output file.
    """
    import numpy as np
    import nilearn.image as niimg
    from   six import string_types

    for arg in list(imgs.keys()):
        if isinstance(imgs[arg], string_types):
            continue

        if np.isscalar(imgs[arg]):
            if arg not in formula:
                raise ValueError("Could not find {} in the formula: {}.".format(arg, formula))

            formula = formula.replace(arg, str(imgs[arg]))
            imgs.pop(arg)

    return niimg.math_img(formula=formula, **imgs)


@ni2file(suffix='_resampled')
def resample(in_file, **kwargs):
    """ Use nilearn.image.resample_img.

    Returns
    -------
    out_file: str
        The absolute path to the output file.
    """
    import nilearn.image as niimg
    return niimg.resample_img(img=in_file, **kwargs)


@ni2file(suffix='_resampled')
def resample_to_img(in_file, target, **kwargs):
    """ Use nilearn.image.resample_to_img.

    Returns
    -------
    out_file: str
        The absolute path to the output file.
    """
    import nilearn.image as niimg
    return niimg.resample_to_img(source_img=in_file, target_img=target, **kwargs)


@ni2file(out_file='concat_img.nii.gz')
def concat_imgs(in_files, out_file=None):
    """ Use nilearn.image.concat_imgs to concat images of up to 4 dimensions.
    Returns
    -------
    out_file: str
        The absolute path to the output file.
    """
    import nilearn.image as niimg
    return niimg.concat_imgs(in_files)


@ni2file(out_file='concat_img.nii.gz')
def concat_3D_imgs(in_files, out_file=None):
    """ Use nilearn.image.concat_imgs to concat 3D volumes into one 4D volume.

    If `in_files` is a list of 3D volumes the return value is the path to one 4D volume.
    Else if `in_files` is a list of 4D volumes the return value is `in_files`.

    Returns
    -------
    out_file: str
        The absolute path to the output file.
    """
    import nilearn.image as niimg

    from   nilearn._utils import check_niimg_3d

    all_3D = True
    for idx, img in enumerate(in_files):
        try:
            _ = check_niimg_3d(img)
        except Exception:
            all_3D = False
            break

    if not all_3D:
        #raise AttributeError('Expected all input images to be 3D volumes, but '
        #                     ' at least the {}th is not.'.format(idx))
        return in_files
    else:
        return niimg.concat_imgs(in_files)


@ni2file(suffix='_mean')
def mean_img(in_file, out_file=None):
    """ Use nilearn.image.mean_img.
    Returns
    -------
    out_file: str
        The absolute path to the output file.
    """
    import nilearn.image as niimg
    return niimg.mean_img(in_file)


@ni2file(suffix='_smooth')
def smooth_img(in_file, fwhm, out_file=None):
    """ Use nilearn.image.smooth_img.
    Returns
    -------
    out_file: str
        The absolute path to the output file.
    """
    import nilearn.image as niimg
    return niimg.smooth_img(in_file, fwhm=fwhm)


@ni2file(suffix='')
def copy_header(in_file, data_file, out_file=None):
    """ Use nilearn.image.new_img_like to copy the header
    from `in_file` to `data_file` and return the result.
    Returns
    -------
    out_file: str
        The absolute path to the output file.
    """
    import nilearn.image as niimg

    img = niimg.load_img(data_file)

    return niimg.new_img_like(in_file, img.get_data(),
                              affine=img.affine, copy_header=True)