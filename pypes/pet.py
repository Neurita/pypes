"""
PET image preprocessing nipype workflows.
"""
import os.path as op

import nipype.pipeline.engine as pe
from   nipype.interfaces.utility import Select
from   nipype.interfaces.base import traits
from   nipype.algorithms.misc import Gunzip

from   .anat import attach_t1_preprocessing
from   .registration import spm_apply_deformations, spm_coregister
from   .petpvc import PETPVC
from   .utils import fsl_merge


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

    #pvc.run()

    return pvc


def spm_pet_preprocessing():
    """ Run the PET-FDG pre-processing workflow against the pet_fdg files in `data_dir`.
    It depends on the t1_preproc_workflow, so if this has not been run, this function
    will run it too.

    It does:
    - SPM12 Coregister PET to T1
    - SPM12 Warp PET to MNI

    Nipype Inputs
    -------------
    gunzip_pet.in_file: traits.File
        The raw NIFTI_GZ PET image file

    coreg_pet.target: traits.File
        Target of the co-registration process, i.e., the anatomical image in native space.

    warp_pet.deformation_file: traits.File
        The DARTEL deformation file for warping the PET to MNI

    tissues_lst.inlist: list of traits.File
        List of tissues files. This node will select only the first 3 images for merging for
        the PET PVC process.

    Returns
    -------
    wf: nipype Workflow
    """
    def flatten_list(list_of_lists):
        if not list_of_lists:
            return []
        if isinstance(list_of_lists[0], list):
            return [l.pop() for l in list_of_lists]
        return list_of_lists

    # fixed parameters of the NUK mMR
    psf_fwhm = (4.3, 4.3, 4.3)

    # coreg pet
    gunzip_pet  = pe.Node(Gunzip(),                             name="gunzip_pet")
    coreg_pet   = pe.Node(spm_coregister(cost_function='mi'),   name="coreg_pet")
    warp_pet    = pe.Node(spm_apply_deformations(),             name="warp_pet")
    pet_pvc     = pe.Node(petpvc_cmd(fwhm_mm=psf_fwhm),         name="pet_pvc")

    anat_to_pet = pe.Node(spm_coregister(cost_function='mi'),   name="anat_to_pet")

    # petpvc
    tissues_lst = pe.Node(Select(),    name="tissues_lst")
    merge_tissu = pe.Node(fsl_merge(), name="merge_tissues")
    tissues_lst.inputs.set(index=[0, 1, 2]) # pick the 3 first items in the list of tissues: GM, WM and CSF
    merge_tissu.inputs.merged_file = "merged_tissues.nii.gz"

    # Create the workflow object
    wf = pe.Workflow(name="pet_preproc")

    # Connect the nodes
    wf.connect([
                (gunzip_pet,  anat_to_pet, [("out_file", "target")]),
                (tissues_lst, anat_to_pet, [(("out", flatten_list), "apply_to_files")]),

                (anat_to_pet, merge_tissu, [("coregistered_files",  "in_files")]),

                (merge_tissu, pet_pvc,     [("merged_file",   "mask_file")]),
                (gunzip_pet,  pet_pvc,     [("out_file",      "in_file")]),

                #(coreg_pet,   warp_pet,    [("coregistered_source", "apply_to_files")]),


                #(gunzip_pet,   anat_to_pet, [("out_file", "target")]),
                #(gunzip_pet,  coreg_pet,   [("out_file",            "source")]),
                #(coreg_pet,   warp_pet,    [("coregistered_source", "apply_to_files")]),
                #(tissues_lst, merge_tissu, [(("out", flatten_list), "in_files")]),

                #(coreg_pet,   pet_pvc,     [("coregistered_source", "in_file")]),
                #(merge_tissu, pet_pvc,     [("merged_file",         "mask_file")]),
              ])

    return wf


def attach_pet_preprocessing(main_wf, data_dir, work_dir=None, output_dir=None):
    """

    Parameters
    ----------
    main_wf
    data_dir
    work_dir
    output_dir

    Returns
    -------

    """
    # Dependency workflows
    main_wf = attach_t1_preprocessing(main_wf=main_wf,
                                      data_dir=data_dir,
                                      work_dir=work_dir,
                                      output_dir=output_dir,)

    input_files = main_wf.get_node("input_files")
    new_segment = main_wf.get_node("t1_preproc.new_segment")
    datasink    = main_wf.get_node("datasink")

    # get the PET preprocessing pipeline
    pet_wf = spm_pet_preprocessing()

    # dataSink output substitutions
    substitutions = [
                     ("wrpet_fdg.nii", "pet_fdg_mni.nii"),
                     ("rpet_fdg.nii",  "pet_fdg_anat.nii"),
                    ]
    datasink.inputs.substitutions.extend(substitutions)

    # Connect the nodes
    main_wf.connect([
                # pet coregistration
                (input_files, pet_wf, [("pet_fdg",                      "gunzip_pet.in_file")]),

                #(new_segment, pet_wf, [("bias_corrected_images",        "coreg_pet.target" )]),
                #(new_segment, pet_wf, [("forward_deformation_field",    "warp_pet.deformation_file")]),

                # pet pvc
                (new_segment, pet_wf, [("native_class_images",          "tissues_lst.inlist")]),
                (new_segment, pet_wf, [("bias_corrected_images",        "anat_to_pet.source")]),

                # datasink
                #(pet_wf,  datasink,   [("coreg_pet.coregistered_files",  "pet.others"),
                #                       ("coreg_pet.coregistered_source", "pet.@anat" )]),
                #(pet_wf,  datasink,   [("warp_pet.normalized_files",     "pet.@mni"  )]),
              ])

    return main_wf
