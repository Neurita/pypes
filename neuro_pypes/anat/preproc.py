# -*- coding: utf-8 -*-
"""
Nipype workflows to process anatomical MRI.
"""
import os.path as op

import nipype.pipeline.engine    as pe
from   nipype.algorithms.misc    import Gunzip
from   nipype.interfaces.utility import IdentityInterface, Function
from   nipype.interfaces.base    import isdefined

from .utils import biasfield_correct, spm_segment

from   ..interfaces.nilearn import math_img, copy_header
from   ..config  import setup_node, check_atlas_file, get_config_setting
from   ..preproc import (spm_apply_deformations,
                         get_bounding_box,)
from   .._utils  import format_pair_list
from   ..utils   import (remove_ext,
                         selectindex,
                         spm_tpm_priors_path,
                         extend_trait_list,
                         get_input_node,
                         get_datasink,
                         get_input_file_name,
                         extension_duplicates)


def spm_anat_preprocessing(wf_name="spm_anat_preproc"):
    """ Run the T1 pre-processing workflow against the anat_hc
    files in `data_dir`.

    It does:
    - N4BiasFieldCorrection
    - SPM12 New Segment
    - SPM12 Warp of MPRAGE to MNI

    [Optional: from config]
    - Atlas file warping to MPRAGE
    - Cortical thickness (SPM+DiReCT)

    Nipype Inputs
    -------------
    anat_input.in_file: traits.File
        Path to the anatomical image.

    anat_input.atlas_file: traits.File
        Path to an atlas file in MNI space to be
        warped to the anatomical space.
        Can also be set through the configuration
        setting `atlas_file`.

    Nipype Outputs
    --------------
    anat_output.anat_mni: traits.File
        The bias-field normalized to MNI anatomical image.

    anat_output.tissues_warped: traits.File
        The tissue segmentation in MNI space from SPM.

    anat_output.tissues_native: traits.File
        The tissue segmentation in native space from SPM.

    anat_output.affine_transform: traits.File
        The affine transformation file.

    anat_output.warp_forward: traits.File
        The forward (anat to MNI) warp field from SPM.

    anat_output.warp_inverse: traits.File
        The inverse (MNI to anat) warp field from SPM.

    anat_output.anat_biascorr: traits.File
        The bias-field corrected anatomical image.

    anat_output.brain_mask: traits.File
        A brain mask file in anatomical space.
        This is calculated by summing up the maps of
        segmented tissues (CSF, WM, GM) and then binarised.

    anat_output.atlas_anat: traits.File
        If `atlas_file` is an existing file in MNI space.
        The atlas file warped to anatomical space,
        if do_atlas and the atlas file is set in configuration.

    anat_output.cortical_thickness: traits.File
        If `anat_preproc.do_cortical_thickness` is True.
        The cortical thickness estimations calculated with the
        SPM+DiReCT method (KellyKapowski).

    anat_output.warped_white_matter: warped_white_matter
        If `anat_preproc.do_cortical_thickness` is True.
        The warped white matter image calculated with the
        SPM+DiReCT method (KellyKapowski).

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
                  "brain_mask",
                 ]

    # check if we have to warp an atlas files too.
    do_atlas, atlas_file = check_atlas_file()
    if do_atlas:
        in_fields  += ["atlas_file"]
        out_fields += ["atlas_anat"]

    # check if we have to do cortical thickness (SPM+DiReCT) method.
    do_cortical_thickness = get_config_setting('anat_preproc.do_cortical_thickness', False)
    if do_cortical_thickness:
        out_fields += ["cortical_thickness",
                       "warped_white_matter",]

    # input node
    anat_input = pe.Node(IdentityInterface(fields=in_fields, mandatory_inputs=True),
                         name="anat_input")

    # atlas registration
    if do_atlas and not isdefined(anat_input.inputs.atlas_file):
        anat_input.inputs.set(atlas_file=atlas_file)

    # T1 preprocessing nodes
    biascor     = setup_node(biasfield_correct(),      name="bias_correction")
    gunzip_anat = setup_node(Gunzip(),                 name="gunzip_anat")
    segment     = setup_node(spm_segment(),            name="new_segment")
    warp_anat   = setup_node(spm_apply_deformations(), name="warp_anat")

    tpm_bbox = setup_node(Function(function=get_bounding_box,
                                   input_names=["in_file"],
                                   output_names=["bbox"]),
                          name="tpm_bbox")
    tpm_bbox.inputs.in_file = spm_tpm_priors_path()

    # calculate brain mask from tissue maps
    tissues = setup_node(IdentityInterface(fields=["gm", "wm", "csf"], mandatory_inputs=True),
                         name="tissues")

    brain_mask = setup_node(Function(function=math_img,
                                     input_names=["formula", "out_file", "gm", "wm", "csf"],
                                     output_names=["out_file"],
                                     imports=['from neuro_pypes.interfaces.nilearn import ni2file']),
                            name='brain_mask')
    brain_mask.inputs.out_file = "tissues_brain_mask.nii.gz"
    brain_mask.inputs.formula  = "np.abs(gm + wm + csf) > 0"

    # output node
    anat_output = pe.Node(IdentityInterface(fields=out_fields), name="anat_output")

    # Connect the nodes
    wf.connect([
                # input to biasfieldcorrection
                (anat_input,   biascor    , [("in_file",      "input_image")]),

                # new segment
                (biascor,      gunzip_anat, [("output_image", "in_file")]),
                (gunzip_anat,  segment,     [("out_file",     "channel_files")]),

                # Normalize12
                (segment,   warp_anat,  [("forward_deformation_field", "deformation_file")]),
                (segment,   warp_anat,  [("bias_corrected_images",     "apply_to_files")]),
                (tpm_bbox,  warp_anat,  [("bbox",                      "write_bounding_box")]),

                # brain mask from tissues
                (segment, tissues,  [(("native_class_images", selectindex, 0), "gm"),
                                     (("native_class_images", selectindex, 1), "wm"),
                                     (("native_class_images", selectindex, 2), "csf"),
                                    ]),

                (tissues, brain_mask, [("gm", "gm"), ("wm", "wm"), ("csf", "csf"),]),

                # output
                (warp_anat, anat_output, [("normalized_files",           "anat_mni")]),
                (segment,   anat_output, [("modulated_class_images",     "tissues_warped"),
                                          ("native_class_images",        "tissues_native"),
                                          ("transformation_mat",         "affine_transform"),
                                          ("forward_deformation_field",  "warp_forward"),
                                          ("inverse_deformation_field",  "warp_inverse"),
                                          ("bias_corrected_images",      "anat_biascorr")]),
                (brain_mask, anat_output, [("out_file",                  "brain_mask")]),
              ])

    # atlas warping nodes
    if do_atlas:
        gunzip_atlas = pe.Node(Gunzip(), name="gunzip_atlas")
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
                    (gunzip_anat,   anat_bbox,    [("out_file",                   "in_file")]),
                    (gunzip_atlas,  warp_atlas,   [("out_file",                   "apply_to_files")]),
                    (segment,       warp_atlas,   [("inverse_deformation_field",  "deformation_file")]),
                    (anat_bbox,     warp_atlas,   [("bbox",                       "write_bounding_box")]),
                    (warp_atlas,    anat_output,  [("normalized_files",           "atlas_anat")]),
                  ])

    # cortical thickness (SPM+DiReCT) method
    if do_cortical_thickness:
        from ..interfaces.ants import KellyKapowski

        segm_img = setup_node(Function(function=math_img,
                                       input_names=["formula", "out_file", "gm", "wm"],
                                       output_names=["out_file"],
                                       imports=['from neuro_pypes.interfaces.nilearn import ni2file']),
                              name='gm-wm_image')
        segm_img.inputs.out_file = "gm_wm.nii.gz"
        segm_img.inputs.formula = '((gm >= 0.5)*2 + (wm > 0.5)*3).astype(np.uint8)'

        # copy the header from the GM tissue image to the result from `gm-wm_image`.
        # this is necessary because the `gm-wm_image` operation sometimes modifies the
        # offset of the image, which will provoke an ANTs exception due to
        # ITK tolerance in ImageToImageFilter
        # https://github.com/stnava/ANTs/issues/74
        cp_hdr = setup_node(Function(function=copy_header,
                                     input_names=["in_file", "data_file"],
                                     output_names=["out_file"],
                                     imports=['from neuro_pypes.interfaces.nilearn import ni2file']),
                            name='copy_header')

        kk = setup_node(KellyKapowski(), name='direct')
        kk.inputs.cortical_thickness  = 'direct_cortical_thickness.nii.gz'
        kk.inputs.warped_white_matter = 'direct_warped_white_matter.nii.gz'

        # connect the cortical thickness (SPM+DiReCT) method
        wf.connect([
                    # create segmentation GM+WM file
                    (tissues, segm_img, [("gm", "gm"),
                                         ("wm", "wm")]),

                    (segm_img, cp_hdr, [("out_file", "data_file")]),
                    (tissues,  cp_hdr, [("gm",       "in_file")]),

                    # kellykapowski
                    (cp_hdr,   kk, [("out_file", "segmentation_image")]),
                    (tissues,  kk, [("gm", "gray_matter_prob_image"),
                                    ("wm", "white_matter_prob_image")]),

                    (kk,    anat_output,  [("cortical_thickness",  "cortical_thickness"),
                                           ("warped_white_matter", "warped_white_matter")]),
                  ])
    return wf


def attach_spm_anat_preprocessing(main_wf, wf_name="spm_anat_preproc"):
    """ Attach the SPM12 anatomical MRI pre-processing workflow to
    the `main_wf`.

    Parameters
    ----------
    main_wf: nipype Workflow

    wf_name: str
        Name of the preprocessing workflow

    Nipype Inputs for `main_wf`
    ---------------------------
    Note: The `main_wf` workflow is expected to have an
    `input_files` and a `datasink` nodes.

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
                    (r"/{anat}_.*corrected_seg8.mat$",    "/{anat}_to_mni_affine.mat"),
                    (r"/m{anat}.*_corrected.nii$",        "/{anat}_biascorrected.nii"),
                    (r"/wm{anat}.*_corrected.nii$",       "/{anat}_mni.nii"),
                    (r"/y_{anat}.*nii$",                  "/{anat}_to_mni_field.nii"),
                    (r"/iy_{anat}.*nii$",                 "/{anat}_to_mni_inv_field.nii"),
                    (r"/mwc1{anat}.*nii$",                "/{anat}_gm_mod_mni.nii"),
                    (r"/mwc2{anat}.*nii$",                "/{anat}_wm_mod_mni.nii"),
                    (r"/mwc3{anat}.*nii$",                "/{anat}_csf_mod_mni.nii"),
                    (r"/mwc4{anat}.*nii$",                "/{anat}_nobrain_mod_mni.nii"),
                    (r"/c1{anat}.*nii$",                  "/{anat}_gm.nii"),
                    (r"/c2{anat}.*nii$",                  "/{anat}_wm.nii"),
                    (r"/c3{anat}.*nii$",                  "/{anat}_csf.nii"),
                    (r"/c4{anat}.*nii$",                  "/{anat}_nobrain.nii"),
                    (r"/c5{anat}.*nii$",                  "/{anat}_nobrain_mask.nii"),
                    (r"/direct_cortical_thickness.nii$",  "/{anat}_gm_cortical_thickness.nii"),
                    (r"/direct_warped_white_matter.nii$", "/{anat}_warped_white_matter.nii"),
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
                     (anat_wf,  datasink, [("anat_output.anat_mni",         "anat.@mni"),
                                           ("anat_output.tissues_warped",   "anat.tissues.warped"),
                                           ("anat_output.tissues_native",   "anat.tissues.native"),
                                           ("anat_output.affine_transform", "anat.transform.@linear"),
                                           ("anat_output.warp_forward",     "anat.transform.@forward"),
                                           ("anat_output.warp_inverse",     "anat.transform.@inverse"),
                                           ("anat_output.anat_biascorr",    "anat.@biascor"),
                                           ("anat_output.brain_mask",       "anat.@brain_mask"),
                                          ]),
                    ])

    # check optional outputs
    if do_atlas:
        main_wf.connect([(anat_wf, datasink, [("anat_output.atlas_anat", "anat.@atlas")]),])

    do_cortical_thickness = get_config_setting('anat_preproc.do_cortical_thickness', False)
    if do_cortical_thickness:
        main_wf.connect([(anat_wf, datasink, [("anat_output.cortical_thickness",  "anat.@cortical_thickness"),
                                              ("anat_output.warped_white_matter", "anat.@warped_white_matter"),
                                             ]),
                        ])

    return main_wf
