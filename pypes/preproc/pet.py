"""
PET image preprocessing nipype function helpers.
"""
import nipype.pipeline.engine    as pe
from   nipype.interfaces         import fsl
from   nipype.interfaces.base    import traits
from   nipype.interfaces.utility import Merge, Split

from   .petpvc import PETPVC
from   ..utils import fsl_merge


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
    split_tissues.inlist: list of existing files
        List of tissue files in anatomical space, the 3 file
        paths must be in this order: GM, WM, CSF

    Nipype.Outputs
    --------------
    merge_tissues.merged_file: existing file
        A 4D volume file with these maps in order: GM, WM, CSF, background

    brain_mask.out_file: existing file
        A mask that is a binarised sum of the tissues file with fslmaths.
        Can be used as brain mask in anatomical space for the PET image.

    Returns
    -------
    wf: nipype Workflow
    """
    # define nodes
    merge_list = pe.Node(Merge(3), name="merge_list")

    split_tissues = pe.Node(Split(splits=[1, 2], squeeze=True), name="split_tissues")

    ## maths for background
    img_bkg = pe.Node(fsl.MultiImageMaths(), name="background")
    img_bkg.inputs.op_string = "-add '%s' -add '%s' -sub 1 -mul -1 -thr 0"
    img_bkg.inputs.out_file  = "tissue_bkg.nii.gz"

    ## maths for brain mask
    brain_mask = pe.Node(fsl.MultiImageMaths(), name="brain_mask")
    brain_mask.inputs.op_string = "-add '%s' -add '%s' -abs -bin"
    brain_mask.inputs.out_file  = "brain_mask.nii.gz"

    merge_tissu = pe.Node(fsl_merge(), name="merge_tissues")
    merge_tissu.inputs.merged_file = "merged_tissues.nii.gz"

    # Create the workflow object
    wf = pe.Workflow(name=wf_name)

    # Connect the nodes
    wf.connect([
                # separate [GM, WM, CSF] into [GM] and [WM, CSF]
                # input [GM] to img_bkg and then [WM, CSF]
                (split_tissues, img_bkg,       [("out1",     "in_file")]),
                (split_tissues, img_bkg,       [("out2",     "operand_files")]),

                (split_tissues, brain_mask,    [("out1",     "in_file")]),
                (split_tissues, brain_mask,    [("out2",     "operand_files")]),

                # create a list of [GM, WM, CSF, BKG]
                (split_tissues, merge_list,    [("out1",     "in1"),
                                                ("out2",     "in2"),
                                               ]),
                (img_bkg,       merge_list,    [("out_file", "in3")]),

                # merge into 4D: [GM, WM, CSF, BKG]
                (merge_list,    merge_tissu,   [("out",       "in_files")]),
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
    mean_value.in_file: existing file
        The image from where to extract the signal values and normalize.

    gm_norm.in_file: existing file
        The image from where to extract the signal values and normalize.
        Usually is the same as mean_value.in_file, but must be explicitly specified anyway.

    mean_value.mask_file: existing file
        The mask to specify which voxels to use to calculate the statistics
        for normalization.

    Nipype Outputs
    --------------
    gm_norm.out_file: existing file

    Returns
    -------
    wf: nipype Workflow
    """
    ## calculate masked stats
    mean_value = pe.Node(fsl.ImageStats(op_string="-M -k %s"), name='mean_value')

    ## normalize
    gm_norm = pe.Node(fsl.BinaryMaths(operation='div'), name='gm_norm')

    # Create the workflow object
    wf = pe.Workflow(name=wf_name)

    wf.connect([
                (mean_value, gm_norm,   [("out_stat",  "operand_value")]),
               ])

    return wf

