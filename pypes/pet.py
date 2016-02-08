# -*- coding: utf-8 -*-
"""
PET image preprocessing nipype workflows.
"""
import os.path as op

import nipype.pipeline.engine    as pe
from   nipype.algorithms.misc    import Gunzip
from   nipype.interfaces.utility import Select, Merge

from   .preproc import spm_apply_deformations, spm_coregister, petpvc_cmd, petpvc_mask, intensity_norm
from   .io      import get_input_file_name
from   .utils   import extend_trait_list, remove_ext, get_input_node, get_datasink
from   ._utils  import flatten_list, format_pair_list


def spm_mrpet_preprocessing(wf_name="spm_mrpet_preproc"):
    """ Run the PET pre-processing workflow against the gunzip_pet.in_file files.
    It depends on the anat_preproc_workflow, so if this has not been run, this function
    will run it too.

    It does:
    - SPM12 Coregister PET to T1
    - SPM12 Warp PET to MNI

    Parameters
    ----------
    wf_name: str
        Name of the workflow.

    Nipype Inputs
    -------------
    gunzip_pet.in_file: traits.File
        The raw NIFTI_GZ PET image file

    coreg_pet.target: traits.File
        Target of the co-registration process, i.e., the anatomical image in native space.

    warp_pet.deformation_file: traits.File
        The DARTEL deformation file for warping the PET to MNI

    tissues.inlist: list of traits.File
        List of tissues files from the New Segment process. At least the first
        3 tissues must be present.

    Nipype outputs
    --------------
    petpvc.out_file: existing file
        The results of the PVC process

    brain_mask.in: existing file
        A brain mask calculated with the tissues file.

    coreg_pet.coregistered_files: list of existing files
        List of coregistered files from coreg_pet.apply_to_files

    coreg_pet.coregistered_source: existing file
        The coregistered PET file

    warp_pet.normalized_files:
        PET images normalized to MNI.
        The result of every internal pre-processing step is normalized to MNI here.

    warp_pet.normalization_parameters: existing files
        Spatial normalization parameters .mat files

    rbvpvc.out_file: existing file
        The result of the PETPVC correction step.

    petpvc_mask.brain_mask.out_file: existing file
        A brain mask calculated through the sum of the 3 tissue files.

    intensity_norm_gm.gm_norm.out_file: existing file
        The output of the Global Mean intensity normalization process.

    Returns
    -------
    wf: nipype Workflow
    """
    # fixed parameters of the NUK mMR
    psf_fwhm = (4.3, 4.3, 4.3)

    # coreg pet
    gunzip_pet  = pe.Node(Gunzip(),                           name="gunzip_pet")
    coreg_pet   = pe.Node(spm_coregister(cost_function="mi"), name="coreg_pet")
    tissues_sel = pe.Node(Select(index=[0, 1, 2]),            name="tissues")
    select_gm   = pe.Node(Select(index=[0]),                  name="select_gm")
    rbvpvc      = pe.Node(petpvc_cmd(fwhm_mm=psf_fwhm,
                                     pvc_method='RBV'),       name="rbvpvc")
    warp_pet    = pe.Node(spm_apply_deformations(),           name="warp_pet")
    merge_lists = pe.Node(Merge(2),                           name='merge_for_warp')

    unzip_mrg = pe.Node(Merge(3),                             name='merge_for_unzip')
    gunzipper = pe.MapNode(Gunzip(),                          name="gunzip", iterfield=['in_file'])

    # workflow to create the mask
    mask_wf = petpvc_mask(wf_name="petpvc_mask")

    # workflow for intensity normalization
    norm_wf = intensity_norm(wf_name="intensity_norm_gm")

    # Create the workflow object
    wf = pe.Workflow(name=wf_name)

    wf.connect([
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
               ])

    return wf


def attach_spm_mrpet_preprocessing(main_wf, wf_name="spm_pet_preproc"):
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
                     (r"/wr{pet}_.*_pvc.nii$",          "/{pet}_mni_pvc.nii.gz"),
                     (r"/wr{pet}_.*_pvc_maths.nii$",    "/{pet}_mni_pvc_norm.nii.gz"),
                     (r"/wbrain_mask.nii",              "/brain_mask_mni.nii"),
                   ]
    regexp_subst = format_pair_list(regexp_subst, pet=pet_fbasename)
    datasink.inputs.regexp_substitutions = extend_trait_list(datasink.inputs.regexp_substitutions,
                                                             regexp_subst)

    # Connect the nodes
    main_wf.connect([
                # pet file input
                (in_files, pet_wf, [("pet",                                  "gunzip_pet.in_file")]),

                # pet to anat registration
                (anat_wf,  pet_wf, [("new_segment.bias_corrected_images",     "coreg_pet.target" )]),
                (anat_wf,  pet_wf, [("new_segment.forward_deformation_field", "warp_pet.deformation_file")]),

                # pet pvc
                (anat_wf,  pet_wf, [("new_segment.native_class_images",       "tissues.inlist")]),

                # datasink
                (pet_wf, datasink, [("coreg_pet.coregistered_files",       "pet.others"),
                                    ("coreg_pet.coregistered_source",      "pet.@anat"),
                                    ("rbvpvc.out_file",                    "pet.@pvc"  ),
                                    ("petpvc_mask.brain_mask.out_file",    "pet.@brain_mask"),
                                    ("intensity_norm_gm.gm_norm.out_file", "pet.@norm"),
                                    ("warp_pet.normalized_files",          "pet.warped2mni"),
                                   ]),
              ])

    return main_wf
