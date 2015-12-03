"""
Nipype workflows to process anatomical MRI.
"""
import os.path as op

import nipype.interfaces.spm     as spm
import nipype.pipeline.engine    as pe
from   nipype.algorithms.misc    import Gunzip
from   nipype.interfaces.ants    import N4BiasFieldCorrection
from   nipype.interfaces.base    import traits
from   nipype.interfaces.io      import DataSink, SelectFiles

from   .preproc import spm_apply_deformations
from   ._utils  import format_pair_list
from   .utils   import (spm_tpm_priors_path,
                        extend_trait_list,
                        find_wf_node,
                        remove_ext)


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

    tissue1 = ((priors_path, 1), 1, (True, True), (True, True))
    tissue2 = ((priors_path, 2), 1, (True, True), (True, True))
    tissue3 = ((priors_path, 3), 2, (True, True), (True, True))
    tissue4 = ((priors_path, 4), 3, (True, True), (True, True))
    tissue5 = ((priors_path, 5), 4, (True, False), (False, False))
    tissue6 = ((priors_path, 6), 2, (False, False), (False, False))
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

    Nipype Inputs
    -------------
    bias_correction.input_image: traits.File
        path to the anatomical image

    Returns
    -------
    wf: nipype Workflow
    """
    # T1 preprocessing nodes
    biascor     = pe.Node(biasfield_correct(),      name="bias_correction")
    gunzip_anat = pe.Node(Gunzip(),                 name="gunzip_anat")
    segment     = pe.Node(spm_segment(),            name="new_segment")
    warp_anat   = pe.Node(spm_apply_deformations(), name="warp_anat")

    # Create the workflow object
    wf = pe.Workflow(name=wf_name)

    # Connect the nodes
    wf.connect([
                # new segment
                (biascor,      gunzip_anat, [("output_image", "in_file"      )]),
                (gunzip_anat,  segment,     [("out_file",     "channel_files")]),

                # Normalize12
                (segment, warp_anat, [("forward_deformation_field", "deformation_file")]),
                (segment, warp_anat, [("bias_corrected_images",     "apply_to_files")]),
              ])
    return wf


def attach_spm_anat_preprocessing(main_wf, data_dir, work_dir=None, output_dir=None, wf_name="spm_anat_preproc"):
    """ Attach the SPM12 anatomical MRI pre-processing workflow to the `main_wf`.

    Parameters
    ----------
    main_wf: nipype Workflow

    data_dir: str

    work_dir: str

    output_dir: str

    wf_name: str
        Name of the anat preprocessing workflow

    Nipype Inputs
    -------------
    Note: The `main_wf` workflow is expected to have an `input_files` and a `datasink` nodes.

    input_files.select.anat: input node

    datasink: nipype Node

    Returns
    -------
    main_wf: nipype Workflow
    """
    input_files = find_wf_node(main_wf, SelectFiles)
    datasink    = find_wf_node(main_wf, DataSink)

    # The workflow box
    t1_wf = spm_anat_preprocessing()

    # The base name of the 'anat' file for the substitutions
    select_node = input_files.get_node('select')
    try:
        anat_fbasename = remove_ext(op.basename(select_node.interface._templates['anat']))
    except:
        raise AttributeError("Could not find a SelectFiles node called 'select' in main workflow.")

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
    datasink.inputs.regexp_substitutions = extend_trait_list(datasink.inputs.regexp_substitutions,
                                                             regexp_subst)

    # input and output anat workflow to main workflow connections
    main_wf.connect([(input_files, t1_wf, [("select.anat",                           "bias_correction.input_image")]),
                     (t1_wf,    datasink, [("warp_anat.normalized_files",            "anat.@mni")],),
                     (t1_wf,    datasink, [("new_segment.modulated_class_images",    "anat.tissues.@warped"),
                                           ("new_segment.native_class_images",       "anat.tissues.@native"),
                                           ("new_segment.transformation_mat",        "anat.transform.@linear"),
                                           ("new_segment.forward_deformation_field", "anat.transform.@forward"),
                                           ("new_segment.inverse_deformation_field", "anat.transform.@inverse"),
                                           ("new_segment.bias_corrected_images",     "anat.@biascor"),
                                          ]),
                    ])

    return main_wf
