"""
PET image preprocessing nipype workflows.
"""
import os.path as op

import nipype.pipeline.engine    as pe
from   nipype.algorithms.misc    import Gunzip
from   nipype.interfaces         import fsl
from   nipype.interfaces.base    import traits
from   nipype.interfaces.utility import Select, Merge, Split
from   nipype.interfaces.io      import DataSink, SelectFiles

from   .anat    import attach_spm_anat_preprocessing
from   .preproc import spm_apply_deformations, spm_coregister, PETPVC
from   .utils   import fsl_merge, extend_trait_list, remove_ext, find_wf_node
from   ._utils  import flatten_list, format_pair_list


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
    split_tissues.inlist: list of existing files
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
    merge_list = pe.Node(Merge(3), name="merge_list")

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
                # input [GM] to img_bkg and then [WM, CSF]
                (split_tissues, img_bkg,       [("out1",     "in_file")]),
                (split_tissues, img_bkg,       [("out2",     "operand_files")]),

                (split_tissues, brain_mask,    [("out1",     "in_file")]),
                (split_tissues, brain_mask,    [("out2",     "operand_files")]),

                # create a list of [GM, WM, CSF, BKG]
                (split_tissues, merge_list,    [("out1",     "in1"),
                                                ("out2",     "in2"),
                                               ]),
                (img_bkg,       merge_list,    [("out_file", "in3")]),

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

    gm_norm.in_file: existing file
        The image from where to extract the signal values and normalize.
        Usually is the same as mean_value.in_file, but must be explicitly specified anyway.

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
    mean_value = pe.Node(fsl.ImageStats(op_string="-M"), name='mean_value')

    ## normalize
    gm_norm = pe.Node(fsl.BinaryMaths(operation='div'), name='gm_norm')

    # Create the workflow object
    wf = pe.Workflow(name=wf_name)

    wf.connect([
                (mean_value, gm_norm,   [("out_stat",  "operand_value")]),
               ])

    return wf


def spm_pet_preprocessing():
    """ Run the PET-FDG pre-processing workflow against the pet_fdg files in `data_dir`.
    It depends on the anat_preproc_workflow, so if this has not been run, this function
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

    warp_pet.normalized_files:
        PET images normalized to MNI

    warp_pet.normalization_parameters: existing files
        Spatial normalization parameters .mat files

    rbvpvc.out_file: existing file

    petpvc_mask.brain_mask.out_file: existing file

    intensity_norm_gm.gm_norm.out_file: existing file


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
    wf = pe.Workflow(name="spm_pet_preproc")

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


def attach_spm_pet_preprocessing(main_wf, data_dir, work_dir=None, output_dir=None):
    """ Attach a PET pre-processing workflow that uses SPM12 to `main_wf`.

    This will also attach the anat preprocessing workflow to `main_wf`. The reason
    for this is that the PET pre-processing steps here make use of anatomical MR
    pre-processing outputs.

    #TODO: a pet pre-processing workflow that does not need anatomical MR
    pre-processing.

    Nipype Inputs
    -------------
    Note: The `main_wf` workflow is expected to have an `input_files` and a `datasink` nodes.

    input_files.select.pet: input node

    datasink: nipype Node

    Parameters
    ----------
    main_wf: nipype Workflow

    data_dir: str

    work_dir: str

    output_dir: str

    Returns
    -------
    main_wf: nipype Workflow
    """
    # Dependency workflows
    main_wf = attach_spm_anat_preprocessing(main_wf=main_wf,
                                            data_dir=data_dir,
                                            work_dir=work_dir,
                                            output_dir=output_dir,
                                            wf_name="spm_anat_preproc")

    anat_wf  = main_wf.get_node("spm_anat_preproc")
    in_files = find_wf_node(main_wf, SelectFiles)
    datasink = find_wf_node(main_wf, DataSink)

    # The base name of the 'pet' file for the substitutions
    select_node = in_files.get_node('select')
    try:
        pet_fbasename = remove_ext(op.basename(select_node.interface._templates['pet']))
    except:
        raise AttributeError("Could not find a SelectFiles node called 'select' in main workflow.")

    # get the PET preprocessing pipeline
    pet_wf = spm_pet_preprocessing()

    # dataSink output substitutions
    regexp_subst = [
                     (r"/r{pet}.nii$",                  "/{pet}_anat.nii"),
                     (r"/r{pet}_.*_pvc.nii.gz$",        "/{pet}_anat_pvc.nii.gz"),
                     (r"/r{pet}_.*_pvc_maths.nii.gz$",  "/{pet}_anat_pvc_norm.nii.gz"),
                     (r"/wr{pet}.nii",                  "/{pet}_mni.nii"),
                     (r"/wr{pet}_.*_pvc.nii.gz$",       "/{pet}_mni_pvc.nii.gz"),
                     (r"/wr{pet}_.*_pvc_maths.nii.gz$", "/{pet}_mni_pvc_norm.nii.gz"),
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
