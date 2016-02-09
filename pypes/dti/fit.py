# -*- coding: utf-8 -*-
"""
Nipype workflows to preprocess diffusion MRI.
"""
import os.path as op
from itertools import product

import nibabel as nib
import numpy   as np
import nipype.pipeline.engine    as pe
from   nipype.interfaces.fsl     import ExtractROI, Eddy, MultiImageMaths
from   nipype.interfaces.io      import DataSink, SelectFiles
from   nipype.interfaces.utility import Function, Select, Split, Merge, IdentityInterface
from   nipype.algorithms.misc    import Gunzip
from   nipype.workflows.dmri.fsl.utils import eddy_rotate_bvecs

from   ..preproc import spm_coregister, spm_apply_deformations
from   ..utils   import find_wf_node, remove_ext, extend_trait_list
from   .._utils  import flatten_list, format_pair_list


def get_bounding_box(in_file):
    """
    Retrieve the bounding box of a volume in millimetres.
    """
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


def write_acquisition_parameters(in_file, epi_factor=128):
    """
    # Comments on the `eddy` tool from FSL FDT.

    A description of the tool:
    http://fsl.fmrib.ox.ac.uk/fsl/fslwiki/Eddy/UsersGuide

    Our problem to run this tool instead of the good-old `eddy_correct` is the `--acqp` argument, an
    acquisitions parameters file.

    A detailed description of the --acpq input file is here:
    http://fsl.fmrib.ox.ac.uk/fsl/fslwiki/eddy/Faq#How_do_I_know_what_to_put_into_my_--acqp_file

    In the following subsections I describe each of the fields needed to check and build the acquisitions parameters file.

    ## Phase Encoding Direction
    Dicom header field: (0018,1312) InPlanePhaseEncodingDirection
    The phase encoding direction is the OPPOSITE of frequency encoding direction:
    - 'COL' = A/P (freqDir L/R),
    - 'ROW' = L/R (freqDir A/P)


    Nifti header field: "descrip.phaseDir:'+'" is for 'COL' in the DICOM InPlanePhaseEncodingDirection value.
    So if you have only one phase encoding oriendation and a '+' in the previous header field,
    the first 3 columns for the `--acqp` parameter file should be:
    0 1 0

    indicating that the scan was acquired with phase-encoding in the anterior-posterior direction.

    For more info:
    http://web.stanford.edu/group/vista/cgi-bin/wiki/index.php/DTI_Preprocessing_User_Manual#Frequency_vs_Phase_Encode_Direction
    https://xwiki.nbirn.org:8443/xwiki/bin/view/Function-BIRN/PhaseEncodeDirectionIssues


    ## Effective Echo Spacing (aka dwell time)

    Effective Echo Spacing (s) = 1/(BandwidthPerPixelPhaseEncode * MatrixSizePhase)

    effective echo spacing = 1 / [(0019,1028) * (0051,100b component #1)] (from the archives)
    https://www.jiscmail.ac.uk/cgi-bin/webadmin?A3=ind1303&L=FSL&E=quoted-printable&P=29358&B=--B_3444933351_15386849&T=text%2Fhtml;%20charset=ISO-8859-1&pending=

    The dwell time is in the nifti header `descrip.dwell`, in seconds (or look at the field `time_units`).
    http://www.mit.edu/~satra/nipype-nightly/interfaces/generated/nipype.interfaces.fsl.epi.html

    More info:
    http://lcni.uoregon.edu/kb-articles/kb-0003

    ## EPI factor

    The EPI factor is not included in the nifti header.
    You can read it using the Grassroots DICOM tool called `gdcmdump`, for example:
    >>> gdcmdump -C IM-0126-0001.dcm | grep 'EPIFactor'
    sFastImaging.lEPIFactor                  = 128

    More info:
    http://dicomlookup.com/default.htm


    # The fourth element of the acquisitions parameter file

    The fourth element in each row is the time (in seconds) between reading the center of the first echo and reading the
    center of the last echo.
    It is the "dwell time" multiplied by "number of PE steps - 1" and it is also the reciprocal of the PE bandwidth/pixel.

    Total readout time (FSL) = (number of echoes - 1) * echo spacing

    Total Readout Time (SPM) = 1/(BandwidthPerPixelPhaseEncode)
    Since the Bandwidth Per Pixel Phase Encode is in Hz, this will give the readout time in seconds


    # The `---index` argument

    The index argument is a text file with a row of numbers. Each number
    indicates what line (starting from 1) in the `acqp` file corresponds to
    each volume in the DTI acquisition.

    # What to do now with the `dcm2nii` files?

    I see two options to calculate the `acqp` lines with these files.

    1. We already have the `dwell` but we don't have the EPI factor.
    We know that the standard in Siemens is 128 and we could stick to that.

    2. Use the `slice_duration * 0.001` which is very near the calculated value.

    # Summary

    So, for example, if we had these acquisition parameters:

    ```
    Phase enc. dir. P >> A
    Echo spacing 0.75 [ms]
    EPI factor 128
    ```

    We should put in the `acqp` file this line:
    0 1 0 0.095
    """
    acqp_file = "diff.acqp"
    index_file = "diff.index"

    image = nib.load(in_file)
    n_directions = image.shape[-1]
    header = image.header
    descrip = dict([item.split("=", 1) for item in header["descrip"][()].split(";")])

    if descrip.get("phaseDir") == "+":
        pe_axis = "0 1 0"
    elif descrip.get("phaseDir") == "-":
        pe_axis = "0 -1 0"
    else:
        raise ValueError("unexpected value for phaseDir: {}".format(descrip.get("phaseDir")))

    # (number of phase-encode steps - 1) * (echo spacing time in milliseconds) * (seconds per millisecond)
    total_readout_time = (epi_factor - 1) * float(descrip["dwell"]) * 1e-3

    with open(acqp_file, "wt") as fout:
        fout.write("{} {}\n".format(pe_axis, total_readout_time))
    with open(index_file, "wt") as fout:
        fout.write("{}\n".format(" ".join(n_directions * ["1"])))

    return op.abspath(acqp_file), op.abspath(index_file)



def fsl_dti_preprocessing(atlas_file, wf_name="fsl_dti_preproc"):
    """ Run the diffusion MRI pre-processing workflow against the diff files in `data_dir`.

    This estimates an affine transform from anat to diff space, applies it to
    the brain mask and an atlas, and performs eddy-current correction.

    Nipype Inputs
    -------------
    dti_input.diff: traits.File
        path to the diffusion MRI image
    dti_input.bval: traits.File
        path to the bvals file
    dti_input.bvec: traits.File
        path to the bvecs file
    dti_input.tissues: traits.File
        paths to the NewSegment c*.nii output files
    dti_input.anat: traits.File
        path to the high-contrast anatomical image
    dti_input.mni_to_anat: traits.File
        path to the warp from MNI space to anat space

    Nipype Workflow Dependencies
    ----------------------------
    This workflow depends on:
    - spm_anat_preproc

    Returns
    -------
    wf: nipype Workflow
    """

    dti_input    = pe.Node(IdentityInterface(
        fields=["diff", "bval", "bvec", "tissues", "anat", "mni_to_anat"],
        mandatory_inputs=True),                                 name="dti_input")
    gunzip_atlas = pe.Node(Gunzip(in_file=atlas_file),          name="gunzip_atlas")
    anat_bbox    = pe.Node(Function(
        input_names=["in_file"],
        output_names=["bbox"],
        function=get_bounding_box),                             name="anat_bbox")
    warp_atlas   = pe.Node(spm_apply_deformations(),            name="warp_atlas")
    write_acqp     = pe.Node(Function(
        input_names=["in_file"],
        output_names=["out_acqp", "out_index"],
        function=write_acquisition_parameters),                 name="write_acqp")
    extract_b0   = pe.Node(ExtractROI(t_min=0, t_size=1),       name="extract_b0")
    gunzip_b0    = pe.Node(Gunzip(),                            name="gunzip_b0")
    coreg_merge  = pe.Node(Merge(2),                            name="coreg_merge")
    coreg_b0     = pe.Node(spm_coregister(cost_function="mi"),  name="coreg_b0")
    brain_sel    = pe.Node(Select(index=[0, 1, 2]),             name="brain_sel")
    coreg_split  = pe.Node(Split(splits=[1, 2, 1], squeeze=True), name="coreg_split")
    brain_merge  = pe.Node(MultiImageMaths(),                   name="brain_merge")
    eddy         = pe.Node(Eddy(),                              name="eddy")
    rot_bvec     = pe.Node(Function(
        input_names=["in_bvec", "eddy_params"],
        output_names=["out_file"],
        function=eddy_rotate_bvecs),                            name="rot_bvec")
    dti_output   = pe.Node(IdentityInterface(
        fields=["diff_corrected", "bvec_rotated", "brain_mask_diff", "atlas_diff"]),
                                                                name="dti_output")

    warp_atlas.inputs.write_interp = 0

    brain_merge.inputs.op_string = "-add '%s' -add '%s' -abs -bin"
    brain_merge.inputs.out_file = "brain_mask_diff_space.nii.gz"

    coreg_b0.inputs.write_interp = 0

    # Create the workflow object
    wf = pe.Workflow(name=wf_name)

    # Connect the nodes
    wf.connect([
                (dti_input,     extract_b0,     [("diff",                   "in_file")]),
                (dti_input,     write_acqp,     [("diff",                   "in_file")]),
                (dti_input,     eddy,           [("diff",                   "in_file")]),
                (dti_input,     eddy,           [("bval",                   "in_bval")]),
                (dti_input,     eddy,           [("bvec",                   "in_bvec")]),
                (dti_input,     rot_bvec,       [("bvec",                   "in_bvec")]),
                (dti_input,     brain_sel,      [("tissues",                "inlist")]),
                (dti_input,     anat_bbox,      [("anat",                   "in_file")]),
                (dti_input,     coreg_b0,       [("anat",                   "source")]),
                (dti_input,     warp_atlas,     [("mni_to_anat",            "deformation_file")]),
                (gunzip_atlas,  warp_atlas,     [("out_file",               "apply_to_files")]),
                (anat_bbox,     warp_atlas,     [("bbox",                   "write_bounding_box")]),
                (write_acqp,    eddy,           [("out_acqp",               "in_acqp"),
                                                 ("out_index",              "in_index")]),
                (extract_b0,    gunzip_b0,      [("roi_file",               "in_file")]),
                (gunzip_b0,     coreg_b0,       [("out_file",               "target")]),
                (brain_sel,     coreg_merge,    [(("out", flatten_list),    "in1")]),
                (warp_atlas,    coreg_merge,    [("normalized_files",       "in2")]),
                (coreg_merge,   coreg_b0,       [("out",                    "apply_to_files")]),
                (coreg_b0,      coreg_split,    [("coregistered_files",     "inlist")]),
                (coreg_split,   brain_merge,    [("out1",                   "in_file")]),
                (coreg_split,   brain_merge,    [("out2",                   "operand_files")]),
                (coreg_split,   dti_output,     [("out3",                   "atlas_diff")]),
                (brain_merge,   eddy,           [("out_file",               "in_mask")]),
                (brain_merge,   dti_output,     [("out_file",               "brain_mask_diff")]),
                (eddy,          dti_output,     [("out_corrected",          "diff_corrected")]),
                (eddy,          rot_bvec,       [("out_parameter",          "eddy_params")]),
                (rot_bvec,      dti_output,     [("out_file",               "bvec_rotated")]),
              ])
    return wf


def attach_fsl_dti_preprocessing(main_wf, wf_name="fsl_dti_preproc", params=None):
    """ Attach the FSL-based diffusion MRI pre-processing workflow to the `main_wf`.

    Parameters
    ----------
    main_wf: nipype Workflow

    wf_name: str
        Name of the preprocessing workflow

    params: dict with parameter values
        atlas_file: str
            Path to the anatomical atlas to be transformed to diffusion MRI space.


    Nipype Inputs for `main_wf`
    ---------------------------
    Note: The `main_wf` workflow is expected to have an `input_files` and a `datasink` nodes.

    input_files.select.diff: input node

    datasink: nipype Node

    Returns
    -------
    main_wf: nipype Workflow
    """
    in_files = find_wf_node(main_wf, SelectFiles)
    datasink = find_wf_node(main_wf, DataSink)
    anat_wf  = main_wf.get_node("spm_anat_preproc")

    atlas_file = config.get("atlas_file", None)
    if atlas_file is None:
        raise ValueError('Expected an existing atlas_file, got {}.'.format(atlas_file))
    if not op.exists(atlas_file):
        raise IOError('Expected an existing atlas_file, got {}.'.format(atlas_file))

    atlas_basename = remove_ext(op.basename(atlas_file))

    # The workflow box
    dti_wf = fsl_dti_preprocessing(atlas_file=atlas_file, wf_name=wf_name)

    regexp_subst = [
                     (r"/rw{atlas}\.nii$", "/{atlas}_diff_space.nii"),
                   ]
    regexp_subst = format_pair_list(regexp_subst, atlas=atlas_basename)
    datasink.inputs.regexp_substitutions = extend_trait_list(datasink.inputs.regexp_substitutions,
                                                             regexp_subst)

    # input and output diffusion MRI workflow to main workflow connections
    main_wf.connect([(in_files, dti_wf,   [("diff",                                  "dti_input.diff"),
                                           ("bval",                                  "dti_input.bval"),
                                           ("bvec",                                  "dti_input.bvec")]),
                     (anat_wf,  dti_wf,   [("new_segment.native_class_images",       "dti_input.tissues"),
                                           ("new_segment.bias_corrected_images",     "dti_input.anat"),
                                           ("new_segment.inverse_deformation_field", "dti_input.mni_to_anat")]),
                     (dti_wf,   datasink, [("dti_output.diff_corrected",             "diff.@eddy_corrected"),
                                           ("dti_output.brain_mask_diff",            "diff.@mask"),
                                           ("dti_output.atlas_diff",                 "diff.@atlas"),
                                           ("dti_output.bvec_rotated",               "diff.@bvec_rotated")]),
                    ])

    return main_wf
