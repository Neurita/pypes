# -*- coding: utf-8 -*-
"""
Function utilities to deal with spatial/coordinates problems in the nifti images.
"""


def get_bounding_box(in_file):
    """ Retrieve the bounding box of a volume in millimetres."""

    # the imports must be inside if you want them to work in a nipype.Function node.
    from itertools import product
    import nibabel as nib
    import numpy   as np

    img = nib.load(in_file)

    # eight corners of the 3-D unit cube [0, 0, 0] .. [1, 1, 1]
    corners = np.array(list(product([0, 1], repeat=3)))
    # scale to the index range of the volume
    corners = corners * (np.array(img.shape[:3]) - 1)
    # apply the affine transform
    corners = img.affine.dot(np.hstack([corners, np.ones((8, 1))]).T).T[:, :3]

    # get the extents
    low_corner  = np.min(corners, axis=0)
    high_corner = np.max(corners, axis=0)

    return [low_corner.tolist(), high_corner.tolist()]