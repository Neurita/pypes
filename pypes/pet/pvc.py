# -*- coding: utf-8 -*-
"""
A PETPVC workflow
"""
import os.path as op

import nipype.pipeline.engine as pe
from   nipype.algorithms.misc import Gunzip
from   nipype.interfaces.utility import Select, IdentityInterface, Function

from  .utils     import (petpvc_cmd,
                         petpvc_mask,
                         intensity_norm)
from   ..config  import setup_node, get_config_setting
from   ..preproc import spm_coregister
from   ..utils   import (get_datasink,
                         extend_trait_list,
                         get_input_node,
                         remove_ext,
                         get_input_file_name,
                         extension_duplicates)

from   .._utils import (flatten_list,
                        format_pair_list)


def petpvc_workflow(wf_name="petpvc"):
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
        List of tissues files from the New Segment process.
        GM

    Nipype outputs
    --------------
    pvc_output.coreg_ref: existing file
        The coregistered reference_file image in PET space.

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

    flat_list = pe.Node(Function(input_names=['list_of_lists'], output_names=['out'],
                                 function=flatten_list),
                        name='flatten_tissue_list')

    # coreg pet
    gunzip_pet  = setup_node(Gunzip(),                           name="gunzip_pet")
    coreg_pet   = setup_node(spm_coregister(cost_function="mi"), name="coreg_pet")

    tissues_sel = setup_node(Select(index=[0, 1, 2]),            name="tissues")
    select_gm   = setup_node(Select(index=[0]),                  name="select_gm")
    rbvpvc      = setup_node(petpvc_cmd(fwhm_mm=psf_fwhm,
                                        pvc_method='RBV'),       name="pvc")

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
                (pet_input,   gunzip_pet,  [("in_file",  "in_file")]),
                (pet_input,   tissues_sel, [("tissues",  "inlist")]),
               ])

    # check how to perform the registration, to decide how to build the pipeline
    anat2pet = get_config_setting('registration.anat2pet', False)

    if anat2pet:
        wf.connect([
                    # inputs
                    (pet_input,   coreg_pet,   [("reference_file", "source")]),

                    # unzip to coregister the reference file (anatomical image) to PET space.
                    (gunzip_pet,  coreg_pet,  [("out_file", "target")]),

                    (tissues_sel, flat_list,  [("out", "list_of_lists")]),
                    (flat_list,   coreg_pet,  [("out", "apply_to_files")]),

                    # the list of tissues to the mask wf and the GM for PET intensity normalization
                    (coreg_pet,   select_gm,  [("coregistered_files", "inlist")]),
                    (coreg_pet,   mask_wf,    [("coregistered_files", "pvcmask_input.tissues")]),

                    # the PET in native space to PVC correction
                    (gunzip_pet,  rbvpvc,     [("out_file", "in_file")]),

                    # the merged file with 4 tissues to PCV correction
                    (mask_wf,     rbvpvc,     [("pvcmask_output.petpvc_mask", "mask_file")]),

                    # normalize voxel values of PET PVCed by demeaning it entirely by
                    # GM PET voxel values
                    (rbvpvc,      norm_wf,    [("out_file", "intnorm_input.source")]),
                    (select_gm,   norm_wf,    [("out",      "intnorm_input.mask")]),

                    # output
                    (coreg_pet,   pvc_output, [("coregistered_source",        "coreg_ref")]),
                    (coreg_pet,   pvc_output, [("coregistered_files",         "coreg_others")]),
                    (rbvpvc,      pvc_output, [("out_file",                   "pvc_out")]),
                    (mask_wf,     pvc_output, [("pvcmask_output.brain_mask",  "brain_mask")]),
                    (mask_wf,     pvc_output, [("pvcmask_output.petpvc_mask", "petpvc_mask")]),
                    (norm_wf,     pvc_output, [("intnorm_output.out_file",    "gm_norm")]),
                   ])
    else: # PET to ANAT
        wf.connect([
                    # inputs
                    (pet_input,   coreg_pet,   [("reference_file",     "target")]),

                    # unzip PET image and set as a source to register it to anatomical space.
                    (gunzip_pet,  coreg_pet,  [("out_file",            "source")]),

                    (tissues_sel, flat_list,  [("out", "list_of_lists")]),
                    (flat_list,   coreg_pet,  [("out", "apply_to_files")]),

                    # the list of tissues to the mask wf and the GM for PET intensity normalization
                    (tissues_sel, select_gm,  [("out", "inlist")]),
                    (flat_list,   mask_wf,    [("out", "pvcmask_input.tissues")]),

                    # the PET in ANAT space to PVC correction
                    (coreg_pet,    rbvpvc,     [("coregistered_source", "in_file")]),

                    # the merged file with 4 tissues to PCV correction
                    (mask_wf,     rbvpvc,     [("pvcmask_output.petpvc_mask", "mask_file")]),

                    # normalize voxel values of PET PVCed by demeaning it entirely by GM PET voxel values
                    (rbvpvc,      norm_wf,    [("out_file", "intnorm_input.source")]),
                    (select_gm,   norm_wf,    [("out",      "intnorm_input.mask")]),

                    # output
                    # TODO: coreg_ref should have a different name in this case
                    (coreg_pet,   pvc_output, [("coregistered_source",         "coreg_ref")]),
                    (coreg_pet,   pvc_output, [("coregistered_files",          "coreg_others")]),
                    (rbvpvc,      pvc_output, [("out_file",                    "pvc_out")]),
                    (mask_wf,     pvc_output, [("pvcmask_output.brain_mask",   "brain_mask")]),
                    (mask_wf,     pvc_output, [("pvcmask_output.petpvc_mask",  "petpvc_mask")]),
                    (norm_wf,     pvc_output, [("intnorm_output.out_file",     "gm_norm")]),
                   ])
    return wf


def no_warp_petpvc(wf_name="petpvc"):
    """ Run the PET pre-processing workflow against the gunzip_pet.in_file files.

    It does:
    - PVC the PET image in PET space

    Parameters
    ----------
    wf_name: str
        Name of the workflow.

    Nipype Inputs
    -------------
    pvc_input.in_file: traits.File
        The PET image file.

    pvc_input.tissues: list of traits.File
        List of tissues files from the in PET space.


    Nipype outputs
    --------------
    pvc_output.pvc_out: existing file
        The output of the PETPVC process.

    pvc_output.petpvc_mask: existing file
        The mask built for the PETPVC.

    pvc_output.brain_mask: existing file
        A brain mask calculated with the tissues file.

    Returns
    -------
    wf: nipype Workflow
    """
    # fixed parameters of the NUK mMR
    psf_fwhm = (4.3, 4.3, 4.3)

    # specify input and output fields
    in_fields  = ["in_file",
                  "tissues"]

    out_fields = ["pvc_out",
                  "petpvc_mask",]

    # input
    pvc_input = setup_node(IdentityInterface(fields=in_fields, mandatory_inputs=True),
                           name="pvc_input")

    # coreg pet
    pvc = setup_node(petpvc_cmd(fwhm_mm=psf_fwhm,
                                pvc_method='RBV'),  name="petpvc")

    # workflow to create the mask
    mask_wf = petpvc_mask(wf_name="petpvc_mask")

    # output
    pvc_output = setup_node(IdentityInterface(fields=out_fields), name="pvc_output")


    # Create the workflow object
    wf = pe.Workflow(name=wf_name)

    wf.connect([
                # the list of tissues to the mask wf and the GM for PET intensity normalization
                (pvc_input, mask_wf, [("tissues", "pvcmask_input.tissues")]),

                # the PET in native space to PVC correction
                (pvc_input,  pvc,    [("in_file", "in_file")]),

                # the merged file with 4 tissues to PCV correction
                (mask_wf, pvc, [("pvcmask_output.petpvc_mask", "mask_file")]),

                # output
                #(coreg_pet,   pvc_output, [("coregistered_source",        "coreg_ref")]),
                #(coreg_pet,   pvc_output, [("coregistered_files",         "coreg_others")]),
                (pvc,         pvc_output, [("out_file",                   "pvc_out")]),
                (mask_wf,     pvc_output, [("pvcmask_output.brain_mask",  "brain_mask")]),
                (mask_wf,     pvc_output, [("pvcmask_output.petpvc_mask", "petpvc_mask")]),
                #(norm_wf,     pvc_output, [("intnorm_output.out_file",    "gm_norm")]),
               ])


    return wf