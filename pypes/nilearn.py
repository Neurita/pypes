# -*- coding: utf-8 -*-
"""
Utilities to use nilearn from nipype
"""


def concat_imgs(in_files, out_file):
    import os.path as op
    import nilearn.image as niimg

    catimg = niimg.concat_imgs(in_files)
    catimg.to_filename(out_file)
    return op.abspath(out_file)