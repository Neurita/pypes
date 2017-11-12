# -*- coding: utf-8 -*-
"""
Helper functions for motion correction of 4D images.
"""
from nipype import Workflow, MapNode
from nipype.interfaces.utility import IdentityInterface
from nipype.interfaces.base import traits
from nipype.interfaces.nipy import SpaceTimeRealigner
from nipype.interfaces.fsl import MCFLIRT, Split, Merge, ApplyWarp

from ..config import setup_node


def nipy_motion_correction(in_file=traits.Undefined):
    """
    Simultaneous motion and slice timing correction algorithm
    If slice_times is not specified, this algorithm performs spatial motion correction
    This interface wraps nipy’s SpaceTimeRealign algorithm [Roche2011] or
    simply the SpatialRealign algorithm when timing info is not provided.

    http://www.mit.edu/~satra/nipype-nightly/interfaces/generated/nipype.interfaces.nipy.preprocess.html#spacetimerealigner

    Returns
    -------
    realign: nipype Interface

    Notes
    -----
    Roche A. A four-dimensional registration algorithm with application to joint correction of motion and slice
    timing in fMRI. IEEE Trans Med Imaging. 2011 Aug;30(8):1546-54. http://dx.doi.org/10.1109/TMI.2011.2131152
    """
    #Run spatial realignment only
    realigner = SpaceTimeRealigner()
    realigner.inputs.in_file = in_file

    #res = realigner.run()
    return realigner


def select_volume(filename, which):
    """Return the middle index of a file
    """
    from nibabel import load
    import numpy as np
    if which.lower() == 'first':
        idx = 0
    elif which.lower() == 'middle':
        idx = int(np.ceil(load(filename).shape[3] / 2))
    else:
        raise Exception('unknown value for volume selection : %s' % which)
    return idx


def fsl_motion_correction(name='realign'):
    """Realign a time series to the middle volume using spline interpolation

    Uses MCFLIRT to realign the time series and ApplyWarp to apply the rigid
    body transformations using spline interpolation (unknown order).

    Nipype Inputs
    -------------
    realign_input.func: exising file
        The path to the fMRI file.

    Nipype Outputs
    --------------
    realign_output.realigned_file: exising file
        The path to the realigned fMRI file.

    realign_output.realign_params: mat file
        .mat file with the affine transformation parameters.

    Example
    -------
    >>> wf = fsl_motion_correction()
    >>> wf.inputs.inputspec.func = 'f3.nii'
    >>> wf.run() # doctest: +SKIP
    """
    wf = Workflow(name=name)

    # input node
    input = setup_node(IdentityInterface(fields=['func']),
                       name='realign_input')
    realigner = setup_node(MCFLIRT(save_mats=True, stats_imgs=True),
                           name='mcflirt')
    splitter  = setup_node(Split(dimension='t'),
                           name='splitter')
    warper    = MapNode(ApplyWarp(interp='spline'),
                        iterfield=['in_file', 'premat'],
                        name='warper')

    joiner    = setup_node(Merge(dimension='t'),
                           name='joiner')

    # output node
    output = setup_node(IdentityInterface(fields=['realigned_file',
                                                  'realign_params',
                                                 ]),
                        name='realign_output')

    wf.connect([
                # input
                (input,     realigner, [("func",                             "in_file"),
                                        (("func", select_volume, 'middle'),  "ref_vol"),
                                       ]),
                # split
                (realigner,  splitter, [("out_file",      "in_file")]),
                (realigner,  warper,   [("mat_file",      "premat"),
                                        ("variance_img",  "ref_file"),
                                       ]),
                # warp
                (splitter,  warper,    [("out_files",      "in_file")]),
                (warper,    joiner,    [("mat_file",       "premat")]),
                # output
                (joiner,    output,    [("merged_file", "realigned_file")]),
                (realigner, output,    [("mat_file",    "realign_params")]),
               ])
    return wf
