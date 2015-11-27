"""
PET image preprocessing nipype workflows.
"""

import nipype.pipeline.engine as pe
from   nipype.interfaces.utility import Select, Merge, Split, IdentityInterface
from   nipype.interfaces.base import traits
from   nipype.algorithms.misc import Gunzip
from   nipype.interfaces import fsl

from   .anat import attach_t1_preprocessing
from   .registration import spm_apply_deformations, spm_coregister
from   .petpvc import PETPVC
from   .utils import fsl_merge, extend_trait_list
from   ._utils import flatten_list


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


def petpvc_mask(wf_name="petpvc_mask"):
    """ A Workflow that returns a 4D merge of 4 volumes for PETPVC: GM, WM, CSF and background.

    Parameters
    ----------
    wf_name: str
        The name of the workflow.

    Nipype.Inputs
    -------------
    tissues.in: list of existing files
        List of tissue files in anatomical space, the 3 file
        paths must be in this order: GM, WM, CSF

    Nipype.Outputs
    --------------
    merge_tissues.merged_file: existing file
        A 4D volume file with these maps in order: GM, WM, CSF, background

    brain_mask.out_file: existing file
        A mask that is a binarised sum of the tissues file with fslmaths.
        Can be used as brain mask in anatomical space for the PET image.

    Returns
    -------
    wf: nipype Workflow
    """
    # define nodes
    merge_list = pe.Node(Merge(2), name="merge_list")

    tissues_lst = pe.Node(IdentityInterface(fields=['in'], mandatory_inputs=True), name="tissues")

    split_tissues = pe.Node(Split(splits=[1, 2], squeeze=True), name="split_tissues")

    ## maths for background
    img_bkg = pe.Node(fsl.MultiImageMaths(), name="background")
    img_bkg.inputs.op_string = "-add '%s' -add '%s' -sub 1 -mul -1 -thr 0"
    img_bkg.inputs.out_file  = "tissue_bkg.nii.gz"

    ## maths for brain mask
    brain_mask = pe.Node(fsl.MultiImageMaths(), name="brain_mask")
    brain_mask.inputs.op_string = "-add '%s' -add '%s' -abs -bin"
    brain_mask.inputs.out_file  = "brain_mask.nii.gz"

    merge_tissu = pe.Node(fsl_merge(), name="merge_tissues")
    merge_tissu.inputs.merged_file = "merged_tissues.nii.gz"

    # Create the workflow object
    wf = pe.Workflow(name=wf_name)

    # Connect the nodes
    wf.connect([
                # separate [GM, WM, CSF] into [GM] and [WM, CSF]
                (tissues_lst,   split_tissues, [("in",       "inlist")]),

                # input to img_bkg [GM] and then [WM, CSF]
                (split_tissues, img_bkg,       [("out1",     "in_file")]),
                (split_tissues, img_bkg,       [("out2",     "operand_files")]),

                (split_tissues, brain_mask,    [("out1",     "in_file")]),
                (split_tissues, brain_mask,    [("out2",     "operand_files")]),

                # create a list of [GM, WM, CSF, BKG]
                (tissues_lst,   merge_list,    [("in",       "in1")]),
                (img_bkg,       merge_list,    [("out_file", "in2")]),

                # merge into 4D: [GM, WM, CSF, BKG]
                (merge_list,    merge_tissu,   [("out",       "in_files")]),
              ])

    return wf


def intensity_norm(wf_name='intensity_norm'):
    """ Workflow that uses a mask against a source from where the mean value will be taken.
    This mean value will be used to demean the whole source and leave it in out_file.

    Parameters
    ----------
    wf_name: str
        The name of the workflow.


    Nipype Inputs
    -------------
    mean_value.in_file: existing file
        The image from where to extract the signal values and normalize.

    mean_value.mask_file: existing file
        The mask to specify which voxels to use to calculate the statistics
        for normalization.

    Nipype Outputs
    --------------
    gm_norm.out_file: existing file

    Returns
    -------
    wf: nipype Workflow
    """
    ## calculate masked stats
    mean_value = pe.Node(fsl.ImageStats(), name='mean_value')
    mean_value.op_string = "-M"

    ## normalize
    gm_norm = pe.Node(fsl.BinaryMaths(operation='div'), name='gm_norm')

    # Create the workflow object
    wf = pe.Workflow(name=wf_name)

    wf.connect([
                # unzip to coregister to anatomical image. 'coreg_pet.target' is an input for this wf.
                (mean_value, gm_norm,   [("in_file",   "in_file")]),
                (mean_value, gm_norm,   [("out_stat",  "operand_value")]),
               ])

    return wf


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

    tissues.inlist: list of traits.File
        List of tissues files from the New Segment process. At least the first
        3 tissues must be present.

    Nipype outputs
    --------------
    petpvc.out_file: existing file
        The results of the PVC process

    pet_norm.out

    brain_mask.in: existing file
        A brain mask calculated with the tissues file.

    coreg_pet.coregistered_files: list of existing files
        List of coregistered files from coreg_pet.apply_to_files

    coreg_pet.coregistered_source: existing file
        The coregistered PET source file

    warp_pet.normalized_source:
        Normalized source files

    warp_pet.normalized_files:
        Normalized other files

    warp_pet.normalization_parameters: existing files
        Spatial normalization parameters .mat files

    rbvpvc.out_file: existing file

    brain_mask.out_file: existing file

    norm_pet.out_file: existing file


    Returns
    -------
    wf: nipype Workflow
    """
    # fixed parameters of the NUK mMR
    psf_fwhm = (4.3, 4.3, 4.3)

    # coreg pet
    gunzip_pet  = pe.Node(Gunzip(),                                 name="gunzip_pet")
    coreg_pet   = pe.Node(spm_coregister(cost_function="mi"),       name="coreg_pet")
    warp_pet    = pe.Node(spm_apply_deformations(),                 name="warp_pet")
    tissues_lst = pe.Node(Select(index=[0, 1, 2]),                  name="tissues")
    select_gm   = pe.Node(Select(index=[0], squeeze=True),          name="select_gm")
    brain_mask  = pe.Node(IdentityInterface(fields=['out_file']),   name="brain_mask")
    norm_pet    = pe.Node(IdentityInterface(fields=['out_file']),   name="norm_pet")
    rbvpvc      = pe.Node(petpvc_cmd(fwhm_mm=psf_fwhm,
                                     pvc_method='RBV'),             name="rbvpvc")

    merge_lists = pe.Node(Merge(4), name='merge_for_warp')

    # workflow to create the mask
    mask_wf = petpvc_mask(wf_name="petpvc_mask")

    # workflow for intensity normalization
    norm_wf = intensity_norm(wf_name="intensity_norm_gm")

    # Create the workflow object
    wf = pe.Workflow(name="pet_preproc")

    wf.add_nodes([tissues_lst])

    wf.connect([
                # unzip to coregister to anatomical image. 'coreg_pet.target' is an input for this wf.
                (gunzip_pet,  coreg_pet,   [("out_file",                  "source")]),

                # the list of tissues to the mask wf and the GM for PET intensity normalization
                (tissues_lst, mask_wf,     [(("out", flatten_list),       "tissues.in")]),
                (tissues_lst, norm_wf,     [(("out", flatten_list),       "mask.in")]),

                # the coregistered PET to PVC correction
                (coreg_pet,   rbvpvc,      [("coregistered_source",       "in_file")]),

                # the merged file with 4 tissues to PCV correction
                (mask_wf,     rbvpvc,      [("merge_tissues.merged_file", "mask_file")]),

                # identityinterface to output the brain_mask
                (mask_wf,     brain_mask,  [("brain_mask.out_file",       "out_file")]),

                # normalize voxel values of PET PVCed by demeaning whole by GM PET voxel values
                (rbvpvc,      norm_wf,     [("out_file",                  "mean_value.in_file")]),
                (select_gm,   norm_wf,     [("out",                       "mean_value.mask_file")]),
                (norm_wf,     norm_pet,    [("gm_norm.out_file",          "out_file")]),

                # warp the PET PVCed to MNI
                (mask_wf,     merge_lists, [("brain_mask.out_file",       "in1")]),
                (tissues_lst, merge_lists, [(("out", flatten_list),       "in2")]),
                (rbvpvc,      merge_lists, [("out_file",                  "in3")]),
                (norm_wf,     merge_lists, [("gm_norm.out_file",          "in4")]),

                (coreg_pet,   warp_pet,    [("coregistered_source",       "source")]),
                (merge_lists, warp_pet,    [("out",                       "apply_to_files")]),
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
    datasink.inputs.substitutions = extend_trait_list(datasink.inputs.substitutions,
                                                      substitutions)

    # Connect the nodes
    main_wf.connect([
                # pet coregistration
                (input_files, pet_wf, [("select.pet_fdg",                "gunzip_pet.in_file")]),

                (new_segment, pet_wf, [("bias_corrected_images",         "coreg_pet.target" )]),
                (new_segment, pet_wf, [("forward_deformation_field",     "warp_pet.deformation_file")]),

                # pet pvc
                (new_segment, pet_wf, [("native_class_images",           "tissues.inlist")]),

                # datasink
                (pet_wf,  datasink,   [("coreg_pet.coregistered_files",  "pet.others"),
                                       ("coreg_pet.coregistered_source", "pet.@anat" ),
                                       ("warp_pet.normalized_source"     "pet.@mni"  ),
                                       ("rbvpvc.out_file",               "pet.@pvc"  ),
                                       ("brain_mask.out_file",           "pet.@brain_mask"),
                                       ("norm_pet.out_file",             "pet.@norm"),
                                      ]),
              ])

    return main_wf
