# -*- coding: utf-8 -*-
"""
PET-MR image preprocessing nipype workflows.
"""
import os.path as op

import nipype.pipeline.engine    as pe
from   nipype.algorithms.misc    import Gunzip
from   nipype.interfaces.utility import Select, Merge, IdentityInterface

from   ..config  import setup_node, check_atlas_file
from   ..preproc import (spm_apply_deformations,
                         spm_coregister,
                         petpvc_cmd,
                         petpvc_mask,
                         intensity_norm)

from   ..utils import (get_datasink,
                       extend_trait_list,
                       get_input_node,
                       remove_ext,
                       get_input_file_name,
                       extension_duplicates)

from   .._utils import (flatten_list,
                        format_pair_list)


def spm_mrpet_preprocessing(wf_name="spm_mrpet_preproc"):
    """ Run the PET pre-processing workflow against the gunzip_pet.in_file files.
    It depends on the anat_preproc_workflow, so if this has not been run, this function
    will run it too.

    It does:
    - PVC the PET image
    - SPM12 Coregister PET to T1
    - SPM12 Warp PET to MNI

    Parameters
    ----------
    wf_name: str
        Name of the workflow.

    Nipype Inputs
    -------------
    pet_input.in_file: traits.File
        The raw NIFTI_GZ PET image file

    pet_input.coreg_target: traits.File
        Target of the co-registration process, i.e., the anatomical image in native space.

    pet_input.warp_field: traits.File
        The deformation file for warping the PET to MNI

    pet_input.tissues: list of traits.File
        List of tissues files from the New Segment process. At least the first
        3 tissues must be present.

    pet_input.atlas_anat: traits.File
        Atlas in anatomical native space.

    Nipype outputs
    --------------
    pet_output.out_file: existing file
        The results of the PVC process

    pet_output.brain_mask: existing file
        A brain mask calculated with the tissues file.

    pet_output.coreg_pet: existing file
        The coregistered PET file

    pet_output.coreg_others: list of existing files
        List of coregistered files from coreg_pet.apply_to_files

    pet_output.mni_pet: existing file
        PET images normalized to MNI.
        The result of every internal pre-processing step is normalized to MNI here.

    pet_output.warp_field: existing files
        Spatial normalization parameters .mat files

    pet_output.gm_norm: existing file
        The output of the Global Mean intensity normalization process.

    pet_output.atlas_pet: existing file
        Atlas image warped to PET space.
        If the `atlas_file` option is an existing file and `normalize_atlas` is True.

    Returns
    -------
    wf: nipype Workflow
    """
    # fixed parameters of the NUK mMR
    psf_fwhm = (4.3, 4.3, 4.3)

    # input
    pet_input = setup_node(IdentityInterface(fields=["in_file", "coreg_target", "warp_field", "tissues"]),
                                             name="pet_input")

    # coreg pet
    gunzip_pet  = setup_node(Gunzip(),                           name="gunzip_pet")
    coreg_pet   = setup_node(spm_coregister(cost_function="mi"), name="coreg_pet")
    tissues_sel = setup_node(Select(index=[0, 1, 2]),            name="tissues")
    select_gm   = setup_node(Select(index=[0]),                  name="select_gm")
    rbvpvc      = setup_node(petpvc_cmd(fwhm_mm=psf_fwhm,
                                        pvc_method='RBV'),       name="rbvpvc")
    warp_pet    = setup_node(spm_apply_deformations(),           name="warp_pet")
    merge_lists = setup_node(Merge(2),                           name='merge_for_warp')

    unzip_mrg = setup_node(Merge(3),                             name='merge_for_unzip')
    gunzipper = pe.MapNode(Gunzip(),                             name="gunzip", iterfield=['in_file'])

    # output
    pet_output = setup_node(IdentityInterface(fields=["out_file",
                                                      "brain_mask",
                                                      "coreg_others",
                                                      "coreg_pet",
                                                      "mni_pet",
                                                      "warp_field",
                                                      "pvc_out",
                                                      "pvc_mask",
                                                      "gm_norm",
                                                      "atlas_pet",]),
                                               name="pet_output")

    # workflow to create the mask
    mask_wf = petpvc_mask(wf_name="petpvc_mask")

    # workflow for intensity normalization
    norm_wf = intensity_norm(wf_name="intensity_norm_gm")

    # Create the workflow object
    wf = pe.Workflow(name=wf_name)

    wf.connect([
                # inputs
                (pet_input,   gunzip_pet,  [("in_file",            "in_file")]),
                (pet_input,   coreg_pet,   [("coreg_target",       "target")]),
                (pet_input,   warp_pet,    [("warp_field",         "deformation_file")]),
                (pet_input,   tissues_sel, [("tissues",            "inlist")]),

                # unzip to coregister to anatomical image. 'coreg_pet.target' is an input for this wf.
                (gunzip_pet,  coreg_pet,  [("out_file",            "source")]),

                # the list of tissues to the mask wf and the GM for PET intensity normalization
                (tissues_sel, select_gm,  [(("out", flatten_list), "inlist")]),
                (tissues_sel, mask_wf,    [(("out", flatten_list), "split_tissues.inlist")]),

                # the coregistered PET to PVC correction
                (coreg_pet,   rbvpvc,     [("coregistered_source", "in_file")]),

                # the merged file with 4 tissues to PCV correction
                (mask_wf,     rbvpvc,     [("merge_tissues.merged_file", "mask_file")]),

                # normalize voxel values of PET PVCed by demeaning whole by GM PET voxel values
                (rbvpvc,      norm_wf,    [("out_file",            "mean_value.in_file")]),
                (rbvpvc,      norm_wf,    [("out_file",            "gm_norm.in_file")]),
                (select_gm,   norm_wf,    [("out",                 "mean_value.mask_file")]),

                # gunzip some files for SPM Normalize12
                (rbvpvc,      unzip_mrg,  [("out_file",            "in1")]),
                (mask_wf,     unzip_mrg,  [("brain_mask.out_file", "in2")]),
                (norm_wf,     unzip_mrg,  [("gm_norm.out_file",    "in3")]),
                (unzip_mrg,   gunzipper,  [("out",                 "in_file")]),

                # warp the PET PVCed to MNI
                (gunzipper,  merge_lists, [("out_file",            "in1")]),
                (coreg_pet,  merge_lists, [("coregistered_source", "in2")]),

                (merge_lists, warp_pet,   [("out",                 "apply_to_files")]),

                # output
                (rbvpvc,      pet_output, [("out_file",            "out_file")]),
                (mask_wf,     pet_output, [("brain_mask.out_file", "brain_mask")]),
                (coreg_pet,   pet_output, [("coregistered_source", "coreg_pet")]),
                (coreg_pet,   pet_output, [("coregistered_files",  "coreg_others")]),
                (norm_wf,     pet_output, [("gm_norm.out_file",    "gm_norm")]),
                (warp_pet,    pet_output, [("normalized_files",    "mni_pet")]),
                (warp_pet,    pet_output, [("deformation_field",   "warp_field")]),
               ])

    # add more nodes if to perform atlas registration
    do_atlas, _ = check_atlas_file()
    if do_atlas:
        coreg_atlas = setup_node(spm_coregister(cost_function="mi"), name="coreg_atlas")

        # set the registration interpolation to nearest neighbour.
        coreg_atlas.inputs.write_interp = 0
        wf.connect([
            (pet_input,   coreg_atlas, [("coreg_target",       "source")]),
            (gunzip_pet,  coreg_atlas, [("out_file",           "target")]),
            (pet_input,   coreg_atlas, [("atlas_anat",         "apply_to_files")]),
            (coreg_atlas, pet_output,  [("coregistered_files", "atlas_pet")]),
        ])

    return wf


def attach_spm_mrpet_preprocessing(main_wf, wf_name="spm_mrpet_preproc"):
    """ Attach a PET pre-processing workflow that uses SPM12 to `main_wf`.
    This workflow needs MRI based

    This will also attach the anat preprocessing workflow to `main_wf`. The reason
    for this is that the PET pre-processing steps here make use of anatomical MR
    pre-processing outputs.

    #TODO: a pet pre-processing workflow that does not make use of an anatomical MR
    pre-processing.

    Nipype Inputs for `main_wf`
    ---------------------------
    Note: The `main_wf` workflow is expected to have an `input_files` and a `datasink` nodes.

    input_files.select.pet: input node

    datasink: nipype Node

    Parameters
    ----------
    main_wf: nipype Workflow

    wf_name: str
        Name of the preprocessing workflow

    Nipype Workflow Dependencies
    ----------------------------
    This workflow depends on:
    - spm_anat_preproc

    Returns
    -------
    main_wf: nipype Workflow
    """
    # Dependency workflows
    anat_wf  = main_wf.get_node("spm_anat_preproc")
    in_files = get_input_node(main_wf)
    datasink = get_datasink  (main_wf)

    # The base name of the 'pet' file for the substitutions
    pet_fbasename = remove_ext(op.basename(get_input_file_name(in_files, 'pet')))

    # get the PET preprocessing pipeline
    pet_wf = spm_mrpet_preprocessing(wf_name=wf_name)

    # dataSink output substitutions
    regexp_subst = [
                     (r"/r{pet}.nii$",                  "/{pet}_anat.nii"),
                     (r"/r{pet}_.*_pvc.nii.gz$",        "/{pet}_anat_pvc.nii.gz"),
                     (r"/r{pet}_.*_pvc_maths.nii.gz$",  "/{pet}_anat_pvc_norm.nii.gz"),
                     (r"/wr{pet}.nii",                  "/{pet}_mni.nii"),
                     (r"/wr{pet}_.*_pvc.nii$",          "/{pet}_mni_pvc.nii"),
                     (r"/wr{pet}_.*_pvc_maths.nii$",    "/{pet}_mni_pvc_norm.nii"),
                     (r"/wbrain_mask.nii",              "/brain_mask_mni.nii"),
                   ]
    regexp_subst = format_pair_list(regexp_subst, pet=pet_fbasename)

    # prepare substitution for atlas_file, if any
    do_atlas, atlas_file = check_atlas_file()
    if do_atlas:
        atlas_basename = remove_ext(op.basename(atlas_file))
        regexp_subst.extend([
                             (r"/[\w]*{atlas}\.nii$", "/{atlas}_fmri_space.nii"),
                            ])
        regexp_subst = format_pair_list(regexp_subst, atlas=atlas_basename, pet=pet_fbasename)

    regexp_subst += extension_duplicates(regexp_subst)
    datasink.inputs.regexp_substitutions = extend_trait_list(datasink.inputs.regexp_substitutions,
                                                             regexp_subst)

    # Connect the nodes
    main_wf.connect([
                # pet file input
                (in_files, pet_wf, [("pet",                                   "pet_input.in_file")]),

                # pet to anat registration
                (anat_wf,  pet_wf, [("new_segment.bias_corrected_images",     "pet_input.coreg_target"),
                                    ("new_segment.forward_deformation_field", "pet_input.warp_field"),
                                    ("new_segment.native_class_images",       "pet_input.tissues"),
                                   ]),

                (pet_wf, datasink, [
                                    ("pet_output.out_file",     "mrpet.@pvc"),
                                    ("pet_output.coreg_others", "mrpet.others"),
                                    ("pet_output.coreg_pet",    "mrpet.@anat"),
                                    ("pet_output.brain_mask",   "mrpet.@brain_mask"),
                                    ("pet_output.gm_norm",      "mrpet.@norm"),
                                    ("pet_output.mni_pet",      "mrpet.warped2mni"),
                                   ]),
              ])

    if do_atlas:
            main_wf.connect([(anat_wf,  pet_wf,   [("anat_output.atlas_warped", "pet_input.atlas_anat")]),
                             (pet_wf,   datasink, [("rest_output.atlas_pet",    "mrpet.@atlas")]),
                            ])

    return main_wf
