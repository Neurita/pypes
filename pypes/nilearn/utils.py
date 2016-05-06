# -*- coding: utf-8 -*-
"""
Utilities to use nilearn from nipype
"""


def resample(in_file, **kwargs):
    """ Use nilearn.image.resample_img."""
    import os.path as op
    from nipype.utils.filemanip import fname_presuffix
    from nilearn.image import resample_img
    res_img = resample_img(img=in_file, **kwargs)

    out_file = op.basename(in_file)
    out_file = fname_presuffix(out_file, suffix='_resampled')
    res_img.to_filename(out_file)

    return op.abspath(out_file)


def resample_to_img(source, target, **kwargs):
    """ Use nilearn.image.resample_to_img."""
    import os.path as op
    from nipype.utils.filemanip import fname_presuffix
    from nilearn.image import resample_to_img

    res_img = resample_to_img(source_img=source, target_img=target, **kwargs)

    out_file = op.basename(source)
    out_file = fname_presuffix(out_file, suffix='_resampled')
    res_img.to_filename(out_file)

    return op.abspath(out_file)


def concat_imgs(in_files, out_file):
    import os.path as op
    import nilearn.image as niimg

    catimg = niimg.concat_imgs(in_files)
    catimg.to_filename(out_file)
    return op.abspath(out_file)
