"""
Nipype workflows to process anatomical MRI.
"""
import nipype.pipeline.engine as pe
from   nipype.algorithms.misc import Gunzip
from   nipype.interfaces.ants import N4BiasFieldCorrection
from   nipype.interfaces.base import traits
import nipype.interfaces.spm as spm

from   .utils import spm_tpm_priors_path, extend_trait_list
from   .registration import spm_apply_deformations


def biasfield_correct(anat_filepath=traits.Undefined):
    """ Inhomogeneity correction. Call N4BiasFieldCorrection throught nipype. """

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
    """

    Parameters
    ----------
    anat_filepath: str

    priors_path: str

    Returns
    -------

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


def spm_t1_preprocessing():
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
    wf = pe.Workflow(name="t1_preproc")

    # Connect the nodes
    wf.connect([
                # new segment
                (biascor,      gunzip_anat, [("output_image", "in_file"      )]),
                (gunzip_anat,  segment,     [("out_file",     "channel_files")]),

                # warp
                (segment, warp_anat,  [("forward_deformation_field", "deformation_file"),
                                       ("bias_corrected_images",     "apply_to_files"),
                                      ]),
              ])
    return wf


def attach_t1_preprocessing(main_wf, data_dir, work_dir=None, output_dir=None):
    """ Attach to `main_wf`

    Parameters
    ----------
    main_wf
    data_dir
    work_dir
    output_dir

    Returns
    -------

    """
    input_files = main_wf.get_node("input_files")
    datasink    = main_wf.get_node("datasink")

    # Dependency workflow
    t1_wf = spm_t1_preprocessing()

    # dataSink output substitutions
    substitutions = [
                     ("manat_hc_corrected.nii", "anat_hc_bc.nii"),
                     ("wanat_hc_bc.nii",        "anat_hc_mni.nii"),
                    ]

    datasink.inputs.substitutions = extend_trait_list(datasink.inputs.substitutions,
                                                      substitutions)

    regexp_subst = [
                     (r"anat_.*corrected_seg8.mat", "anat_to_mni_affine.mat"),
                     (r"/y_anat.*nii$",   "/anat_to_mni_field.nii"),
                     (r"/iy_anat.*nii$",  "/anat_to_mni_inv_field.nii"),
                     (r"/mwc1anat.*nii$", "/ant_mni_gm_mod.nii"),
                     (r"/mwc2anat.*nii$", "/anat_mni_wm_mod.nii"),
                     (r"/mwc3anat.*nii$", "/anat_mni_csf_mod.nii"),
                     (r"/mwc4anat.*nii$", "/anat_mni_nobrain_mod.nii"),
                     (r"/c1anat.*nii$",   "/anat_gm.nii"),
                     (r"/c2anat.*nii$",   "/anat_wm.nii"),
                     (r"/c3anat.*nii$",   "/anat_csf.nii"),
                     (r"/c4anat.*nii$",   "/anat_nobrain.nii"),
                     (r"/c5anat.*nii$",   "/anat_nobrain_mask.nii")
                    ]
    datasink.inputs.regexp_substitutions = extend_trait_list(datasink.inputs.regexp_substitutions,
                                                             regexp_subst)

    main_wf.connect([(input_files, t1_wf, [("select.anat_hc",  "bias_correction.input_image")]),

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
