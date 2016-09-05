# -*- coding: utf-8 -*-
"""
A PETPVC workflow
"""
import os.path as op

import nipype.pipeline.engine    as pe
from   nipype.algorithms.misc    import Gunzip
from   nipype.interfaces.utility import Select, Merge, IdentityInterface

from  .utils     import (petpvc_cmd,
                         petpvc_mask,
                         intensity_norm)
from   ..config  import setup_node
from   ..preproc import spm_coregister
from   ..utils   import (get_datasink,
                         extend_trait_list,
                         get_input_node,
                         remove_ext,
                         get_input_file_name,
                         extension_duplicates)

from   .._utils import (flatten_list,
                        format_pair_list)


def petpvc_workflow(wf_name="pet_pvc"):
    """ Run the PET pre-processing workflow against the gunzip_pet.in_file files.
    It coregisters the reference_file and tissues to PET space, then applies PVC and grey matter normalization.

    It does:
    - SPM12 Coregister T1 and tisues to PET
    - PVC the PET image in PET space

    Parameters
    ----------
    wf_name: str
        Name of the workflow.

    Nipype Inputs
    -------------
    pvc_input.in_file: traits.File
        The raw NIFTI_GZ PET image file

    pvc_input.reference_file: traits.File
        The anatomical image in its native space. For registration reference.

    pvc_input.tissues: list of traits.File
        List of tissues files from the New Segment process. At least the first
        3 tissues must be present.

    Nipype outputs
    --------------
    pvc_output.coreg_ref: existing file
        The coregistered reference image to PET space.

    pvc_output.coreg_others: list of existing files
        List of coregistered files from coreg_pet.apply_to_files

    pvc_output.pvc_out: existing file
        The output of the PETPVC process.

    pvc_output.petpvc_mask: existing file
        The mask built for the PETPVC.

    pvc_output.brain_mask: existing file
        A brain mask calculated with the tissues file.

    pvc_output.gm_norm: existing file
        The output of the grey matter intensity normalization process.
        This is the last step in the PET signal correction.

    Returns
    -------
    wf: nipype Workflow
    """
    # fixed parameters of the NUK mMR
    psf_fwhm = (4.3, 4.3, 4.3)

    # specify input and output fields
    in_fields  = ["in_file",
                  "reference_file",
                  "tissues",]

    out_fields = ["coreg_ref",
                  "coreg_others",
                  "pvc_out",
                  "petpvc_mask",
                  "brain_mask",
                  "gm_norm",]

    # input
    pet_input = setup_node(IdentityInterface(fields=in_fields, mandatory_inputs=True),
                           name="pvc_input")

    # coreg pet
    gunzip_pet  = setup_node(Gunzip(),                           name="gunzip_pet")
    coreg_pet   = setup_node(spm_coregister(cost_function="mi"), name="coreg_pet")

    tissues_sel = setup_node(Select(index=[0, 1, 2]),            name="tissues")
    select_gm   = setup_node(Select(index=[0]),                  name="select_gm")
    rbvpvc      = setup_node(petpvc_cmd(fwhm_mm=psf_fwhm,
                                        pvc_method='RBV'),       name="rbvpvc")

    unzip_mrg = setup_node(Merge(3),                             name='merge_for_unzip')
    gunzipper = pe.MapNode(Gunzip(),                             name="gunzip", iterfield=['in_file'])

    # output
    pvc_output = setup_node(IdentityInterface(fields=out_fields), name="pvc_output")

    # workflow to create the mask
    mask_wf = petpvc_mask(wf_name="petpvc_mask")

    # workflow for intensity normalization
    norm_wf = intensity_norm(wf_name="intensity_norm_gm")

    # Create the workflow object
    wf = pe.Workflow(name=wf_name)

    wf.connect([
                # inputs
                (pet_input,   gunzip_pet,  [("in_file",            "in_file")]),
                (pet_input,   coreg_pet,   [("reference_file",     "source")]),
                (pet_input,   tissues_sel, [("tissues",            "inlist")]),

                # unzip to coregister the reference file (anatomical image) to PET space.
                (gunzip_pet,  coreg_pet,  [("out_file",            "target")]),
                (tissues_sel, coreg_pet,  [(("out", flatten_list), "apply_to_files")]),

                # the list of tissues to the mask wf and the GM for PET intensity normalization
                (coreg_pet,   select_gm,  [("coregistered_files", "inlist")]),
                (coreg_pet,   mask_wf,    [("coregistered_files", "pvcmask_input.tissues")]),

                # the PET in native space to PVC correction
                (gunzip_pet,  rbvpvc,     [("out_file", "in_file")]),

                # the merged file with 4 tissues to PCV correction
                (mask_wf,     rbvpvc,     [("pvcmask_output.petpvc_mask", "mask_file")]),

                # normalize voxel values of PET PVCed by demeaning it entirely by GM PET voxel values
                (rbvpvc,      norm_wf,    [("out_file", "intnorm_input.source")]),
                (select_gm,   norm_wf,    [("out",      "intnorm_input.mask")]),

                # gunzip some files for SPM Normalize12
                (rbvpvc,      unzip_mrg,  [("out_file",                    "in1")]),
                (mask_wf,     unzip_mrg,  [("pvcmask_output.brain_mask",   "in2")]),
                (norm_wf,     unzip_mrg,  [("intnorm_output.out_file",     "in3")]),
                (unzip_mrg,   gunzipper,  [("out",                         "in_file")]),

                # output
                (coreg_pet,   pvc_output, [("coregistered_source",         "coreg_ref")]),
                (coreg_pet,   pvc_output, [("coregistered_files",          "coreg_others")]),
                (rbvpvc,      pvc_output, [("out_file",                    "pvc_out")]),
                (mask_wf,     pvc_output, [("pvcmask_output.brain_mask",   "brain_mask")]),
                (mask_wf,     pvc_output, [("pvcmask_output.petpvc_mask",  "petpvc_mask")]),
                (norm_wf,     pvc_output, [("intnorm_output.out_file",     "gm_norm")]),
               ])

    return wf


def attach_petpvc_workflow(main_wf, wf_name="spm_petpvc"):
    """ Attach a PETPVC workflow.

    This will also attach the anat preprocessing workflow to `main_wf`. The reason
    for this is that the PET pre-processing steps here make use of anatomical MR
    pre-processing outputs.

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
    anat_fbasename = remove_ext(op.basename(get_input_file_name(in_files, 'anat')))
    pet_fbasename  = remove_ext(op.basename(get_input_file_name(in_files, 'pet')))

    # get the PET preprocessing pipeline
    pet_wf = petpvc_workflow(wf_name=wf_name)

    # dataSink output substitutions
    # dataSink output substitutions
    regexp_subst = [
                     (r"/{pet}_.*_pvc.nii.gz$",       "/{pet}_pvc.nii.gz"),
                     (r"/{pet}_.*_pvc_maths.nii.gz$", "/{pet}_pvc_norm.nii.gz"),
                     (r"/rm{anat}_corrected.nii$",    "/{anat}_{pet}.nii"),
                     (r"/rc1{anat}_corrected.nii$",   "/gm_{pet}.nii"),
                     (r"/rc2{anat}_corrected.nii$",   "/wm_{pet}.nii"),
                     (r"/rc3{anat}_corrected.nii$",   "/csf_{pet}.nii"),
                   ]
    regexp_subst = format_pair_list(regexp_subst, pet=pet_fbasename, anat=anat_fbasename)
    regexp_subst += extension_duplicates(regexp_subst)
    datasink.inputs.regexp_substitutions = extend_trait_list(datasink.inputs.regexp_substitutions,
                                                             regexp_subst)

    # Connect the nodes
    main_wf.connect([
                     # pet file input
                     (in_files, pet_wf, [("pet", "pvc_input.in_file")]),

                     # pet to anat registration
                     (anat_wf,  pet_wf, [("new_segment.bias_corrected_images", "pet_input.reference_file"),
                                         ("new_segment.native_class_images",   "pet_input.tissues"),
                                        ]),

                     (pet_wf, datasink, [
                                         ("pvc_output.coreg_others", "pet.tissues"),
                                         ("pvc_output.coreg_ref",    "pet.@anat"),
                                         ("pvc_output.pvc_out",      "mrpet.@pvc"),
                                         ("pvc_output.petpvc_mask",  "mrpet.@petpvc_mask"),
                                         ("pvc_output.brain_mask",   "pet.@brain_mask"),
                                         ("pvc_output.gm_norm",      "pet.@gm_norm"),
                                        ]),
                     ])

    return main_wf