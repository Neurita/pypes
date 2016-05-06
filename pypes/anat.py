# -*- coding: utf-8 -*-
"""
Nipype workflows to process anatomical MRI.
"""
import os.path as op

import nipype.interfaces.spm     as spm
import nipype.pipeline.engine    as pe
from   nipype.algorithms.misc    import Gunzip
from   nipype.interfaces.ants    import N4BiasFieldCorrection
from   nipype.interfaces.base    import traits
from   nipype.interfaces.utility import IdentityInterface, Function

from   .config  import setup_node, check_atlas_file
from   .preproc import (spm_apply_deformations,
                        get_bounding_box,)
from   ._utils  import format_pair_list
from   .utils   import (remove_ext,
                        spm_tpm_priors_path,
                        extend_trait_list,
                        get_input_node,
                        get_datasink,
                        get_input_file_name,
                        extension_duplicates)


def biasfield_correct(anat_filepath=traits.Undefined):
    """ Inhomogeneity correction.
    ANTS N4BiasFieldCorrection interface.

    Parameters
    ----------
    anat_filepath: str
        Path to the anatomical file path

    Returns
    -------
    seg: N4BiasFieldCorrection interface
    """
    n4 = N4BiasFieldCorrection()
    n4.inputs.dimension = 3
    n4.inputs.bspline_fitting_distance = 300
    n4.inputs.shrink_factor = 3
    n4.inputs.n_iterations = [50, 50, 30, 20]
    n4.inputs.convergence_threshold = 1e-6
    #n4.inputs.bspline_order = 5
    n4.inputs.save_bias = True

    n4.inputs.input_image = anat_filepath

    return n4


def spm_segment(anat_filepath=traits.Undefined, priors_path=None):
    """ SPM12 New Segment interface.

    Parameters
    ----------
    anat_filepath: str
        Path to the anatomical file path

    priors_path: str
        Path to the tissue probability maps file

    Returns
    -------
    seg: NewSegment interface
    """
    if priors_path is None:
        priors_path = spm_tpm_priors_path()

    seg = spm.NewSegment()

    tissue1 = ((priors_path, 1), 1, (True,  True),   (True,  True))
    tissue2 = ((priors_path, 2), 1, (True,  True),   (True,  True))
    tissue3 = ((priors_path, 3), 2, (True,  True),   (True,  True))
    tissue4 = ((priors_path, 4), 3, (True,  True),   (True,  True))
    tissue5 = ((priors_path, 5), 4, (True,  False),  (False, False))
    tissue6 = ((priors_path, 6), 2, (False, False),  (False, False))
    seg.inputs.tissues = [tissue1, tissue2, tissue3, tissue4, tissue5, tissue6]
    seg.inputs.channel_info = (0.0001, 60, (True, True))
    #seg.inputs.warping_regularization = [0, 0.001, 0.5, 0.05, 0.2]
    seg.inputs.write_deformation_fields = [True, True]

    seg.inputs.channel_files = anat_filepath

    #seg.run()
    return seg


def spm_anat_preprocessing(wf_name="spm_anat_preproc"):
    """ Run the T1 pre-processing workflow against the anat_hc files in `data_dir`.

    It does:
    - N4BiasFieldCorrection
    - SPM12 New Segment
    - SPM12 Warp of MPRAGE to MNI

    [Optional: from config]


    Nipype Inputs
    -------------
    anat_input.in_file: traits.File
        path to the anatomical image

    Nipype Outputs
    --------------
    anat_output.anat_mni: traits.File
        The bias-field normalized to MNI anatomical image.

    anat_output.tissues_warped: traits.File
        The tissue segmentation in MNI space from SPM.

    anat_output.tissues_native: traits.File
        The tissue segmentation in native space from SPM

    anat_output.affine_transform: traits.File
        The affine transformation file.

    anat_output.warp_forward: traits.File
        The forward (anat to MNI) warp field from SPM.

    anat_output.warp_inverse: traits.File
        The inverse (MNI to anat) warp field from SPM.

    anat_output.anat_biascorr: traits.File
        The bias-field corrected anatomical image

    anat_output.atlas_anat: traits.File
        The atlas file warped to anatomical space,
        if do_atlas and the atlas file is set in configuration.

    Returns
    -------
    wf: nipype Workflow
    """
    # Create the workflow object
    wf = pe.Workflow(name=wf_name)

    # specify input and output fields
    in_fields  = ["in_file"]
    out_fields = ["anat_mni",
                  "tissues_warped",
                  "tissues_native",
                  "affine_transform",
                  "warp_forward",
                  "warp_inverse",
                  "anat_biascorr",
                 ]

    do_atlas, atlas_file = check_atlas_file()
    if do_atlas:
        in_fields  += ["atlas_file"]
        out_fields += ["atlas_anat"]

    # input node
    anat_input = setup_node(IdentityInterface(fields=in_fields, mandatory_inputs=True),
                            name="anat_input")
    # atlas registration
    if do_atlas:
        anat_input.inputs.set(atlas_file=atlas_file)

    # T1 preprocessing nodes
    biascor     = setup_node(biasfield_correct(),      name="bias_correction")
    gunzip_anat = setup_node(Gunzip(),                 name="gunzip_anat")
    segment     = setup_node(spm_segment(),            name="new_segment")
    warp_anat   = setup_node(spm_apply_deformations(), name="warp_anat")

    # output node
    anat_output = setup_node(IdentityInterface(fields=out_fields),
                             name="anat_output")

    # Connect the nodes
    wf.connect([
                # input
                (anat_input,   biascor    , [("in_file",      "input_image")]),
                # new segment
                (biascor,      gunzip_anat, [("output_image", "in_file")]),
                (gunzip_anat,  segment,     [("out_file",     "channel_files")]),

                # Normalize12
                (segment,   warp_anat,  [("forward_deformation_field", "deformation_file")]),
                (segment,   warp_anat,  [("bias_corrected_images",     "apply_to_files")]),

                # output
                (warp_anat, anat_output, [("normalized_files",           "anat_mni")]),
                (segment,   anat_output, [("modulated_class_images",     "tissues_warped"),
                                          ("native_class_images",        "tissues_native"),
                                          ("transformation_mat",         "affine_transform"),
                                          ("forward_deformation_field",  "warp_forward"),
                                          ("inverse_deformation_field",  "warp_inverse"),
                                          ("bias_corrected_images",      "anat_biascorr")]),
              ])

    # atlas warping nodes
    if do_atlas:
        gunzip_atlas = setup_node(Gunzip(),                 name="gunzip_atlas")
        warp_atlas   = setup_node(spm_apply_deformations(), name="warp_atlas")
        anat_bbox    = setup_node(Function(function=get_bounding_box,
                                           input_names=["in_file"],
                                           output_names=["bbox"]),
                                  name="anat_bbox")

        # set the warping interpolation to nearest neighbour.
        warp_atlas.inputs.write_interp = 0

        # connect the atlas registration nodes
        wf.connect([
                    (anat_input,    gunzip_atlas, [("atlas_file",                 "in_file")]),
                    (anat_input,    anat_bbox,    [("in_file",                    "in_file")]),
                    (gunzip_atlas,  warp_atlas,   [("out_file",                   "apply_to_files")]),
                    (segment,       warp_atlas,   [("inverse_deformation_field",  "deformation_file")]),
                    (anat_bbox,     warp_atlas,   [("bbox",                       "write_bounding_box")]),
                    (warp_atlas,    anat_output,  [("normalized_files",           "atlas_anat")]),
                  ])
    return wf


def attach_spm_anat_preprocessing(main_wf, wf_name="spm_anat_preproc"):
    """ Attach the SPM12 anatomical MRI pre-processing workflow to the `main_wf`.

    Parameters
    ----------
    main_wf: nipype Workflow

    wf_name: str
        Name of the preprocessing workflow

    Nipype Inputs for `main_wf`
    ---------------------------
    Note: The `main_wf` workflow is expected to have an `input_files` and a `datasink` nodes.

    input_files.anat: input node

    datasink: nipype Node

    Returns
    -------
    main_wf: nipype Workflow
    """
    in_files = get_input_node(main_wf)
    datasink = get_datasink  (main_wf)

    # The workflow box
    anat_wf = spm_anat_preprocessing(wf_name=wf_name)

    # The base name of the 'anat' file for the substitutions
    anat_fbasename = remove_ext(op.basename(get_input_file_name(in_files, 'anat')))

    # dataSink output substitutions
    regexp_subst = [
                    (r"/{anat}_.*corrected_seg8.mat$", "/{anat}_to_mni_affine.mat"),
                    (r"/m{anat}.*_corrected.nii$",     "/{anat}_biascorrected.nii"),
                    (r"/wm{anat}.*_corrected.nii$",    "/{anat}_mni.nii"),
                    (r"/y_{anat}.*nii$",               "/{anat}_to_mni_field.nii"),
                    (r"/iy_{anat}.*nii$",              "/{anat}_to_mni_inv_field.nii"),
                    (r"/mwc1{anat}.*nii$",             "/{anat}_gm_mod_mni.nii"),
                    (r"/mwc2{anat}.*nii$",             "/{anat}_wm_mod_mni.nii"),
                    (r"/mwc3{anat}.*nii$",             "/{anat}_csf_mod_mni.nii"),
                    (r"/mwc4{anat}.*nii$",             "/{anat}_nobrain_mod_mni.nii"),
                    (r"/c1{anat}.*nii$",               "/{anat}_gm.nii"),
                    (r"/c2{anat}.*nii$",               "/{anat}_wm.nii"),
                    (r"/c3{anat}.*nii$",               "/{anat}_csf.nii"),
                    (r"/c4{anat}.*nii$",               "/{anat}_nobrain.nii"),
                    (r"/c5{anat}.*nii$",               "/{anat}_nobrain_mask.nii"),
                   ]
    regexp_subst = format_pair_list(regexp_subst, anat=anat_fbasename)

    # prepare substitution for atlas_file, if any
    do_atlas, atlas_file = check_atlas_file()
    if do_atlas:
        atlas_basename = remove_ext(op.basename(atlas_file))
        regexp_subst.extend([
                             (r"/w{atlas}\.nii$", "/{atlas}_anat_space.nii"),
                            ])
        regexp_subst = format_pair_list(regexp_subst, atlas=atlas_basename)

    # add nii.gz patterns
    regexp_subst += extension_duplicates(regexp_subst)
    datasink.inputs.regexp_substitutions = extend_trait_list(datasink.inputs.regexp_substitutions,
                                                             regexp_subst)

    main_wf.connect([(in_files, anat_wf,  [("anat",                         "anat_input.in_file")]),
                     (anat_wf,  datasink, [("anat_output.anat_mni",         "anat.@mni")],),
                     (anat_wf,  datasink, [("anat_output.tissues_warped",   "anat.tissues.@warped"),
                                           ("anat_output.tissues_native",   "anat.tissues.@native"),
                                           ("anat_output.affine_transform", "anat.transform.@linear"),
                                           ("anat_output.warp_forward",     "anat.transform.@forward"),
                                           ("anat_output.warp_inverse",     "anat.transform.@inverse"),
                                           ("anat_output.anat_biascorr",    "anat.@biascor"),
                                          ]),
                    ])

    # check optional outputs
    if do_atlas:
        main_wf.connect([(anat_wf, datasink, [("anat_output.atlas_anat", "anat.@atlas")]),])

    return main_wf
