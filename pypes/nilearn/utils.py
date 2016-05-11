# -*- coding: utf-8 -*-
"""
Utilities to use nilearn from nipype
"""
from functools import wraps


def ni2file(**presuffixes):
    """ Pick the nibabel image container output from `f` and stores it in a file.
    To know the path where it has to save the file it looks into argument
    values in this order:
    - `out_file` kwarg in the call to `f`,
    - `out_file` kwarg in the `ni2file` decorator,
    - the first argument in the function call.
    In the last case a presuffix must be defined in the decorator to avoid overwriting
    an existing file.
    """
    def nifti_out(f):
        @wraps(f)
        def wrapped(*args, **kwargs):
            import os.path as op
            from nipype.utils.filemanip import fname_presuffix

            res_img = f(*args, **kwargs)

            out_file = kwargs.get('out_file', presuffixes.pop('out_file', None))
            if out_file is None:
                out_file = args[1]
                if not op.exists(out_file):
                    raise IOError('Expected an existing file to use as reference for'
                                  ' the output file name, got {}.'.format(out_file))

                if not presuffixes:
                    raise IOError('The file {} already exists, please add a presuffix to the'
                                  'decorator.'.format(out_file))

            out_file = fname_presuffix(out_file, **presuffixes)
            res_img.to_filename(out_file)

            return op.abspath(out_file)

        return wrapped
    return nifti_out


@ni2file(out_file='nilearn_maths.nii.gz')
def math_img(formula, out_file='', **imgs):
    """ Use nilearn.image.math_img.

    Returns
    -------
    out_file: str
        The absolute path to the output file.
    """
    from nilearn.image import math_img
    return math_img(formula=formula, **imgs)


@ni2file(suffix='_resampled')
def resample(in_file, **kwargs):
    """ Use nilearn.image.resample_img.

    Returns
    -------
    out_file: str
        The absolute path to the output file.
    """
    from nilearn.image import resample_img
    return resample_img(img=in_file, **kwargs)


@ni2file(suffix='_resampled')
def resample_to_img(source, target, **kwargs):
    """ Use nilearn.image.resample_to_img.

    Returns
    -------
    out_file: str
        The absolute path to the output file.
    """
    from nilearn.image import resample_to_img
    return resample_to_img(source_img=source, target_img=target, **kwargs)


@ni2file(out_file='concat_img.nii.gz')
def concat_imgs(in_files, out_file=''):
    """
    Returns
    -------
    out_file: str
        The absolute path to the output file.
    """
    import nilearn.image as niimg
    return niimg.concat_imgs(in_files)
