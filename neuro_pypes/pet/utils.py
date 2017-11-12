# -*- coding: utf-8 -*-
"""
PET image preprocessing utilities nipype function helpers.
"""
from   nipype.interfaces.base import traits
from   nipype.interfaces.utility import Merge, Function, IdentityInterface
from   nipype.pipeline import Workflow

from   ..interfaces.nilearn import math_img, concat_imgs, resample_to_img
from   ..config  import setup_node
from   ..preproc import PETPVC
from   ..utils   import selectindex, rename


#TODO: add a becquerel/ml normalization function node
# hint: img*slope + intercept


def petpvc_cmd(in_file=traits.Undefined, mask_file=traits.Undefined, out_file=traits.Undefined,
               pvc_method='RBV', fwhm_mm=(4, 4, 4)):
    """ Return a nipype interface to PETPVC.

    Parameters
    ----------
    in_file: str
        Path to the PET image file.

    out_file: str
        Path to the result output file

    mask_file: str
        Path to the mask image file with tissue maps.

    pvc_method: str
        Keyword of the method to run.
        Run `petpvc --help` for choices.

    fwhm_mm: iterable of floats
        Iterable with 3 ints or floats to define the full-width at half maximum in mm along
        each axis of the point-spread function (PSF) of the scanner.

    Returns
    -------
    pet_pvc: PETPVC(nipype.interfaces.base.CommandLine)
    """

    pvc = PETPVC()
    pvc.inputs.in_file   = in_file
    pvc.inputs.out_file  = out_file
    pvc.inputs.mask_file = mask_file
    pvc.inputs.pvc       = pvc_method
    pvc.inputs.fwhm_x    = fwhm_mm[0]
    pvc.inputs.fwhm_y    = fwhm_mm[1]
    pvc.inputs.fwhm_z    = fwhm_mm[2]

    return pvc


def petpvc_mask(wf_name="petpvc_mask"):
    """ A Workflow that returns a 4D merge of 4 volumes for PETPVC: GM, WM, CSF and background.

    Parameters
    ----------
    wf_name: str
        The name of the workflow.

    Nipype.Inputs
    -------------
    pvcmask_input.tissues: list of existing files
        List of tissue files in anatomical space, the 3 file
        paths must be in this order: GM, WM, CSF

    Nipype.Outputs
    --------------
    pvcmask_output.petpvc_mask: existing file
        A 4D volume file with these maps in order: GM, WM, CSF, background

    pvcmask_output.brain_mask: existing file
        A mask that is a binarised sum of the tissues file with fslmaths.
        Can be used as brain mask in anatomical space for the PET image.

    Returns
    -------
    wf: nipype Workflow
    """
    # define nodes
    # specify input and output fields
    in_fields  = ["tissues"]

    out_fields = ["petpvc_mask",
                  "brain_mask",]

    # input
    pvcmask_input = setup_node(IdentityInterface(fields=in_fields, mandatory_inputs=True),
                               name="pvcmask_input")

    tissues = setup_node(IdentityInterface(fields=["gm", "wm", "csf"], mandatory_inputs=True),
                         name="tissues")

    merge_list = setup_node(Merge(4), name="merge_list")

    ## maths for background
    img_bkg = setup_node(Function(function=math_img,
                                  input_names=["formula", "out_file", "gm", "wm", "csf"],
                                  output_names=["out_file"],
                                  imports=['from neuro_pypes.interfaces.nilearn import ni2file']),
                          name='background')
    img_bkg.inputs.out_file = "tissue_bkg.nii.gz"
    img_bkg.inputs.formula  = "np.maximum((-((gm + wm + csf) - 1)), 0)"

    ## maths for brain mask
    brain_mask = setup_node(Function(function=math_img,
                                     input_names=["formula", "out_file", "gm", "wm", "csf"],
                                     output_names=["out_file"],
                                     imports=['from neuro_pypes.interfaces.nilearn import ni2file']),
                            name='brain_mask')
    brain_mask.inputs.out_file = "tissues_brain_mask.nii.gz"
    brain_mask.inputs.formula  = "np.abs(gm + wm + csf) > 0"

    ## concat the tissues images and the background for PETPVC
    merge_tissues = setup_node(Function(function=concat_imgs,
                                        input_names=["in_files"],
                                        output_names=["out_file"],
                                        imports=['from neuro_pypes.interfaces.nilearn import ni2file']),
                               name='merge_tissues')
    merge_tissues.inputs.out_file = "petpvc_mask.nii.gz"

    # output
    pvcmask_output = setup_node(IdentityInterface(fields=out_fields), name="pvcmask_output")

    # Create the workflow object
    wf = Workflow(name=wf_name)

    # Connect the nodes
    wf.connect([
                # separate [GM, WM, CSF] into [GM] and [WM, CSF]
                (pvcmask_input, tissues, [(("tissues", selectindex, 0), "gm"),
                                          (("tissues", selectindex, 1), "wm"),
                                          (("tissues", selectindex, 2), "csf"),
                                         ]),

                (tissues, img_bkg,    [("gm", "gm" ), ("wm", "wm" ), ("csf", "csf"),]),
                (tissues, brain_mask, [("gm", "gm" ), ("wm", "wm" ), ("csf", "csf"),]),
                (tissues, merge_list, [("gm", "in1"), ("wm", "in2"), ("csf", "in3"),]),

                # create a list of [GM, WM, CSF, BKG]
                (img_bkg, merge_list, [("out_file", "in4")]),

                # merge into 4D: [GM, WM, CSF, BKG]
                (merge_list, merge_tissues, [("out", "in_files")]),

                # output
                (merge_tissues, pvcmask_output, [("out_file", "petpvc_mask")]),
                (brain_mask,    pvcmask_output, [("out_file", "brain_mask")]),
              ])

    return wf


def intensity_norm(wf_name='intensity_norm'):
    """ Workflow that uses a mask against a source from where the mean value will be taken.
    This mean value will be used to demean the whole source and leave it in out_file.

    Parameters
    ----------
    wf_name: str
        The name of the workflow.

    Nipype Inputs
    -------------
    intnorm_input.source: existing file
        The image from where to extract the signal values and normalize.

    intnorm_input.mask: existing file
        The mask to specify which voxels to use to calculate the statistics
        for normalization.

    Nipype Outputs
    --------------
    intnorm_output.out_file: existing file

    Returns
    -------
    wf: nipype Workflow
    """
    # specify input and output fields
    in_fields  = ["source",
                  "mask"]

    out_fields = ["out_file"]

    # input
    intnorm_input = setup_node(IdentityInterface(fields=in_fields, mandatory_inputs=True),
                               name="intnorm_input")

    # fix the affine matrix (it's necessary for some cases)
    resample = setup_node(Function(function=resample_to_img,
                                   input_names=["in_file", "target", "interpolation"],
                                   output_names=["out_file"],
                                   imports=['from neuro_pypes.interfaces.nilearn import ni2file']),
                          name="resample_mask")
    resample.inputs.interpolation = "nearest"

    # calculate masked mean value
    mean_val = setup_node(Function(function=math_img,
                                   input_names=["formula", "img", "mask"],
                                   output_names=["out_value"],
                                   imports=['from neuro_pypes.interfaces.nilearn import ni2file']),
                          name='mean_value')
    mean_val.inputs.formula = "np.mean(np.nonzero(img[mask > 0]))"

    # normalize
    norm_img = setup_node(Function(function=math_img,
                                   input_names=["formula", "out_file", "img", "val"],
                                   output_names=["out_file"],
                                   imports=['from neuro_pypes.interfaces.nilearn import ni2file']),
                          name='norm_img')
    norm_img.inputs.formula = "img / val"

    # output
    intnorm_output = setup_node(IdentityInterface(fields=out_fields),
                                name="intnorm_output")

    # Create the workflow object
    wf = Workflow(name=wf_name)

    wf.connect([
                # resample
                (intnorm_input, resample,  [("source",    "target"),
                                            ("mask",      "in_file")]),

                # normalize
                (intnorm_input, mean_val,  [("source",    "img" )]),
                (resample,      mean_val,  [("out_file",  "mask")]),

                (intnorm_input, norm_img,  [("source",    "img"),
                                            (("source", rename, "_intnormed"), "out_file"),
                                           ]),

                (mean_val,      norm_img,  [("out_value", "val")]),
                (norm_img, intnorm_output, [("out_file",  "out_file")]),
               ])

    return wf

