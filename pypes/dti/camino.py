# -*- coding: utf-8 -*-
"""
Nipype workflows to use Camino for tractography.
"""
import nipype.pipeline.engine    as pe
import nipype.algorithms.misc    as misc
from   nipype.interfaces.utility import IdentityInterface
from   nipype.interfaces.camino  import (Image2Voxel, FSL2Scheme, DTIFit, Track,
                                         Conmat, ComputeFractionalAnisotropy, AnalyzeHeader)

from   ..utils import (get_datasink,
                       get_input_node,
                       get_data_dims,
                       get_vox_dims,
                       get_affine,
                       )


def camino_tractography(wf_name="camino_tract", fa_tract_stat='mean'):
    """ Run the diffusion MRI pre-processing workflow against the diff files in `data_dir`.

    Nipype Inputs
    -------------
    tract_input.diff: traits.File
        path to the diffusion MRI image
    tract_input.bval: traits.File
        path to the bvals file
    tract_input.bvec: traits.File
        path to the bvecs file
    tract_input.mask: traits.File
        path to the brain mask file
    tract_input.atlas: traits.File
        path to the atlas file

    Returns
    -------
    wf: nipype Workflow
    """

    tract_input  = pe.Node(IdentityInterface(
        fields=["diff", "bvec", "bval", "mask", "atlas"],
        mandatory_inputs=True),                                 name="tract_input")
    img2vox_diff = pe.Node(Image2Voxel(out_type="float"),       name="img2vox_diff")
    img2vox_mask = pe.Node(Image2Voxel(out_type="short"),       name="img2vox_mask")
    fsl2scheme   = pe.Node(FSL2Scheme(),                        name="fsl2scheme")
    dtifit       = pe.Node(DTIFit(),                            name="dtifit")
    fa           = pe.Node(ComputeFractionalAnisotropy(),       name="fa")

    analyzehdr_fa = pe.Node(interface=AnalyzeHeader(), name="analyzeheader_fa")
    analyzehdr_fa.inputs.datatype = "double"
    fa2nii = pe.Node(interface=misc.CreateNifti(), name='fa2nii')

    track        = pe.Node(Track(
        inputmodel="dt",
        out_file="tracts.Bfloat"),                              name="track")
    conmat       = pe.Node(Conmat(output_root="conmat_"),       name="conmat")
    tract_output = pe.Node(IdentityInterface(
        fields=["tensor", "tracks", "connectivity", "mean_fa"]),
                                                                name="tract_output")

    conmat.inputs.tract_stat = fa_tract_stat

    # Create the workflow object
    wf = pe.Workflow(name=wf_name)

    # Connect the nodes
    wf.connect([
                (tract_input,   img2vox_diff,     [("diff",                  "in_file"     )]),
                (tract_input,   fsl2scheme,       [("bvec",                  "bvec_file"   ),
                                                   ("bval",                  "bval_file"   )]),
                (tract_input,   track,            [("atlas",                 "seed_file"   )]),
                (tract_input,   conmat,           [("atlas",                 "target_file" )]),
                (tract_input,   img2vox_mask,     [("mask",                  "in_file"     )]),
                (img2vox_diff,  dtifit,           [("voxel_order",           "in_file"     )]),
                (img2vox_mask,  dtifit,           [("voxel_order",           "bgmask"      )]),
                (fsl2scheme,    dtifit,           [("scheme",                "scheme_file" )]),
                (fsl2scheme,    fa,               [("scheme",                "scheme_file" )]),
                (dtifit,        fa,               [("tensor_fitted",         "in_file"     )]),
                (dtifit,        tract_output,     [("tensor_fitted",         "tensor"      )]),
                (dtifit,        track,            [("tensor_fitted",         "in_file"     )]),
                (track,         conmat,           [("tracked",               "in_file"     )]),
                (track,         tract_output,     [("tracked",               "tracks"      )]),
                (fa,            analyzehdr_fa,    [("fa",                    "in_file"     )]),
                (tract_input,   analyzehdr_fa,    [(('diff', get_vox_dims),  "voxel_dims"  ),
                                                   (('diff', get_data_dims), "data_dims"   )]),
                (analyzehdr_fa, fa2nii,           [("header",                "header_file" )]),
                (tract_input,   fa2nii,           [(("diff", get_affine),    "affine"      )]),
                (fa,            fa2nii,           [("fa",                    "data_file"   )]),
                (fa2nii,        conmat,           [("nifti_file",            "scalar_file" )]),
                (conmat,        tract_output,     [("conmat_sc",             "connectivity"),
                                                   ("conmat_ts",             "mean_fa"     )]),
              ])
    return wf


def attach_camino_tractography(main_wf, wf_name="camino_tract", params=None):
    """ Attach the Camino-based tractography workflow to the `main_wf`.

    Parameters
    ----------
    main_wf: nipype Workflow

    atlas_file: str
        Path to the anatomical atlas.

    wf_name: str
        Name of the preprocessing workflow

    Nipype Inputs for `main_wf`
    ---------------------------
    Note: The `main_wf` workflow is expected to have an `input_files` and a `datasink` nodes.

    input_files.select.diff: input node

    datasink: nipype Node

    Nipype Workflow Dependencies
    ----------------------------
    This workflow depends on:
    - spm_anat_preproc
    - fsl_dti_preproc

    Returns
    -------
    main_wf: nipype Workflow
    """
    in_files = get_input_node(main_wf)
    datasink = get_datasink  (main_wf)
    dti_wf   = main_wf.get_node("fsl_dti_preproc")

    # The workflow box
    tract_wf = camino_tractography(wf_name=wf_name)

    # input and output diffusion MRI workflow to main workflow connections
    main_wf.connect([(in_files, tract_wf, [("bval",                       "tract_input.bval")]),
                     (dti_wf,   tract_wf, [("dti_output.diff_corrected",  "tract_input.diff"),
                                           ("dti_output.bvec_rotated",    "tract_input.bvec"),
                                           ("dti_output.brain_mask_diff", "tract_input.mask"),
                                           ("dti_output.atlas_diff",      "tract_input.atlas")]),
                     (tract_wf, datasink, [("tract_output.tensor",        "tract.@tensor"),
                                           ("tract_output.tracks",        "tract.@tracks"),
                                           ("tract_output.connectivity",  "tract.@connectivity"),
                                           ("tract_output.mean_fa",       "tract.@mean_fa")])
                    ])

    return main_wf
