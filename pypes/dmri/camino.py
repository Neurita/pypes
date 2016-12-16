# -*- coding: utf-8 -*-
"""
Nipype workflows to use Camino for tractography.
"""
import nipype.pipeline.engine    as pe
import nipype.algorithms.misc    as misc
from   nipype.interfaces.utility import IdentityInterface
from   nipype.interfaces.camino  import (Image2Voxel, FSL2Scheme, DTIFit, Track,
                                         Conmat, ComputeFractionalAnisotropy, AnalyzeHeader)

from   ..config  import setup_node, check_atlas_file
from   ..utils import (get_datasink,
                       get_interface_node,
                       get_input_node,
                       get_data_dims,
                       get_vox_dims,
                       get_affine,
                       )


def camino_tractography(wf_name="camino_tract"):
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

    Nipypte Outputs
    ---------------
    tract_output.tensor
        The result of fitting the tensor model to the whole image.

    tract_output.tracks
        The tractography result.

    tract_output.connectivity
        The atlas ROIxROI structural connectivity matrix.

    tract_output.mean_fa
        The atlas ROIxROI structural connectivity matrix with average FA values.

    tract_output.fa
        The voxelwise fractional anisotropy image.

    Returns
    -------
    wf: nipype Workflow
    """
    in_fields  = ["diff", "bvec", "bval", "mask", "atlas"]
    out_fields = ["tensor", "tracks", "connectivity", "mean_fa", "fa"]

    tract_input  = pe.Node(IdentityInterface(fields=in_fields,
                                             mandatory_inputs=True),
                           name="tract_input")

    img2vox_diff = setup_node(Image2Voxel(out_type="float"), name="img2vox_diff")
    img2vox_mask = setup_node(Image2Voxel(out_type="short"), name="img2vox_mask")
    fsl2scheme   = setup_node(FSL2Scheme(),                  name="fsl2scheme")
    dtifit       = setup_node(DTIFit(),                      name="dtifit")
    fa           = setup_node(ComputeFractionalAnisotropy(), name="fa")

    analyzehdr_fa = setup_node(interface=AnalyzeHeader(), name="analyzeheader_fa")
    analyzehdr_fa.inputs.datatype = "double"

    fa2nii = setup_node(interface=misc.CreateNifti(), name='fa2nii')

    track  = setup_node(Track(inputmodel="dt", out_file="tracts.Bfloat"), name="track")
    conmat = setup_node(Conmat(output_root="conmat_"), name="conmat")

    tract_output = pe.Node(IdentityInterface(fields=out_fields),
                           name="tract_output")

    # Create the workflow object
    wf = pe.Workflow(name=wf_name)

    # Connect the nodes
    wf.connect([
                # convert data to camino format
                (tract_input,   img2vox_diff,     [("diff",                  "in_file"     )]),
                (tract_input,   img2vox_mask,     [("mask",                  "in_file"     )]),

                # convert bvec and bval to camino scheme
                (tract_input,   fsl2scheme,       [("bvec",                  "bvec_file"   ),
                                                   ("bval",                  "bval_file"   )]),

                # dtifit
                (img2vox_diff,  dtifit,           [("voxel_order",           "in_file"     )]),
                (img2vox_mask,  dtifit,           [("voxel_order",           "bgmask"      )]),
                (fsl2scheme,    dtifit,           [("scheme",                "scheme_file" )]),

                # calculate FA
                (fsl2scheme,    fa,               [("scheme",                "scheme_file" )]),
                (dtifit,        fa,               [("tensor_fitted",         "in_file"     )]),

                # tractography
                (tract_input,   track,            [("atlas",                 "seed_file"   )]),
                (dtifit,        track,            [("tensor_fitted",         "in_file"     )]),

                # convert FA data to NifTI
                (fa,            analyzehdr_fa,    [("fa",                    "in_file"     )]),
                (tract_input,   analyzehdr_fa,    [(('diff', get_vox_dims),  "voxel_dims"  ),
                                                   (('diff', get_data_dims), "data_dims"   )]),

                (tract_input,   fa2nii,           [(("diff", get_affine),    "affine"      )]),
                (analyzehdr_fa, fa2nii,           [("header",                "header_file" )]),
                (fa,            fa2nii,           [("fa",                    "data_file"   )]),

                # connectivity matrix
                (tract_input,   conmat,           [("atlas",                 "target_file" )]),
                (track,         conmat,           [("tracked",               "in_file"     )]),
                (fa2nii,        conmat,           [("nifti_file",            "scalar_file" )]),

                # output
                (fa2nii,        tract_output,     [("nifti_file",            "fa"          )]),
                (dtifit,        tract_output,     [("tensor_fitted",         "tensor"      )]),
                (track,         tract_output,     [("tracked",               "tracks"      )]),
                (conmat,        tract_output,     [("conmat_sc",             "connectivity"),
                                                   ("conmat_ts",             "mean_fa"     )]),
              ])
    return wf


def attach_camino_tractography(main_wf, wf_name="camino_tract"):
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
    - spm_fsl_dti_preprocessing

    Returns
    -------
    main_wf: nipype Workflow
    """
    in_files = get_input_node(main_wf)
    datasink = get_datasink  (main_wf)
    dti_coreg_output = get_interface_node(main_wf, 'dti_co_output')
    dti_artif_output = get_interface_node(main_wf, 'dti_art_output')

    # The workflow box
    tract_wf = camino_tractography(wf_name=wf_name)

    # input and output diffusion MRI workflow to main workflow connections
    main_wf.connect([(in_files,         tract_wf, [("bval",            "tract_input.bval")]),
                     (dti_coreg_output, tract_wf, [("brain_mask_diff", "tract_input.mask")]),

                     (dti_artif_output, tract_wf, [("eddy_corr_file",  "tract_input.diff"),
                                                   ("bvec_rotated",    "tract_input.bvec"),
                                                  ]),

                     # output
                     (tract_wf, datasink, [("tract_output.tensor",       "tract.@tensor"),
                                           ("tract_output.tracks",       "tract.@tracks"),
                                           ("tract_output.connectivity", "tract.@connectivity"),
                                           ("tract_output.mean_fa",      "tract.@mean_fa"),
                                           ("tract_output.fa",           "tract.@fa"),
                                           ])
                    ])

    # pass the atlas if it's the case
    do_atlas, _ = check_atlas_file()
    if do_atlas:
        main_wf.connect([(dti_coreg_output, tract_wf, [("atlas_diff", "tract_input.atlas")])])

    return main_wf
