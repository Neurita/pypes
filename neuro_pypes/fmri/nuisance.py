# -*- coding: utf-8 -*-
"""
Resting-state fMRI specific nuisance correction filtering workflow.
"""

import nipype.pipeline.engine       as pe
from   nipype.algorithms.rapidart   import ArtifactDetect
from   nipype.interfaces.utility    import Function, IdentityInterface, Merge
from   nipype.algorithms.confounds  import TSNR
from   nipype.interfaces            import fsl

from   ..config  import setup_node, _get_params_for
from   ..utils   import selectindex, rename
from   ..preproc import motion_regressors, extract_noise_components, create_regressors


def rapidart_fmri_artifact_detection():
    art = ArtifactDetect()
    art.inputs.use_differences      = [True, False]
    art.inputs.use_norm             = True
    art.inputs.zintensity_threshold = 2
    art.inputs.norm_threshold       = 1
    art.inputs.mask_type            = 'file'
    art.inputs.parameter_source     = 'NiPy'
    return art


def rest_noise_filter_wf(wf_name='rest_noise_removal'):
    """ Create a resting-state fMRI noise removal node.

    Nipype Inputs
    -------------
    rest_noise_input.in_file

    rest_noise_input.brain_mask

    rest_noise_input.wm_mask

    rest_noise_input.csf_mask

    rest_noise_input.motion_params
        Nipy motion parameters.

    Nipype Outputs
    --------------
    rest_noise_output.tsnr_file
        A SNR estimation volume file for QA purposes.

    rest_noise_output.motion_corrected
        The fMRI motion corrected image.

    rest_noise_output.nuis_corrected
        The resulting nuisance corrected image.
        This will be the same as 'motion_corrected' if compcor
        is disabled.

    rest_noise_output.motion_regressors
        Motion regressors file.

    rest_noise_output.compcor_regressors
        CompCor regressors file.

    rest_noise_output.art_displacement_files
        One image file containing the voxel-displacement timeseries.

    rest_noise_output.art_intensity_files
        One file containing the global intensity values determined
        from the brainmask.

    rest_noise_output.art_norm_files
        One file containing the composite norm.

    rest_noise_output.art_outlier_files
         One file containing a list of 0-based indices corresponding
         to outlier volumes.

    rest_noise_output.art_plot_files
        One image file containing the detected outliers.

    rest_noise_output.art_statistic_files
        One file containing information about the different types of
        artifacts and if design info is provided then details of
        stimulus correlated motion and a listing or artifacts by
        event type.

    Returns
    -------
    rm_nuisance_wf: nipype Workflow
    """

    # Create the workflow object
    wf = pe.Workflow(name=wf_name)

    in_fields = ["in_file",
                 "brain_mask",
                 "wm_mask",
                 "csf_mask",
                 "motion_params",]

    out_fields = ["tsnr_file",
                  "motion_corrected",
                  "nuis_corrected",
                  "motion_regressors",
                  "compcor_regressors",
                  "gsr_regressors",
                  "art_displacement_files",
                  "art_intensity_files",
                  "art_norm_files",
                  "art_outlier_files",
                  "art_plot_files",
                  "art_statistic_files",
                 ]

    # input identities
    rest_noise_input = setup_node(IdentityInterface(fields=in_fields,
                                                    mandatory_inputs=True),
                                  name="rest_noise_input")

    # get the settings for filters
    filters = _get_params_for('rest_filter')

    # Compute TSNR on realigned data regressing polynomial up to order 2
    tsnr = setup_node(TSNR(regress_poly=2), name='tsnr')

    # Use :class:`nipype.algorithms.rapidart` to determine which of the
    # images in the functional series are outliers based on deviations in
    # intensity or movement.
    art = setup_node(rapidart_fmri_artifact_detection(), name="detect_artifacts")

    # Compute motion regressors
    motion_regs = setup_node(Function(input_names=['motion_params',
                                                   'order',
                                                   'derivatives',
                                                  ],
                                      output_names=['out_files'],
                                      function=motion_regressors,),
                             name='motion_regressors')

    # Create a filter to remove motion and art confounds
    motart_pars = setup_node(Function(input_names=['motion_params',
                                                   'comp_norm',
                                                   'outliers',
                                                   'detrend_poly'],
                                      output_names=['out_files'],
                                      function=create_regressors),
                             name='motart_parameters')

    motion_filter = setup_node(fsl.GLM(out_f_name='F_mcart.nii.gz',
                                       out_pf_name='pF_mcart.nii.gz',
                                       demean=True),
                               name='motion_filter')

    # Noise confound regressors
    compcor_pars = setup_node(Function(input_names=['realigned_file',
                                                    'mask_file',
                                                    'num_components',
                                                    'extra_regressors'],
                                         output_names=['out_files'],
                                         function=extract_noise_components,),
                              name='compcor_pars')
    #compcor_pars = setup_node(ACompCor(), name='compcor_pars')
    #compcor_pars.inputs.components_file = 'noise_components.txt'

    compcor_filter = setup_node(fsl.GLM(out_f_name='F.nii.gz',
                                        out_pf_name='pF.nii.gz',
                                        demean=True),
                                name='compcor_filter')

    # Global signal regression
    gsr_pars = setup_node(Function(input_names=['realigned_file',
                                                'mask_file',
                                                'num_components',
                                                'extra_regressors'],
                                    output_names=['out_files'],
                                    function=extract_noise_components, ),
                            name='gsr_pars')

    gsr_filter = setup_node(fsl.GLM(out_f_name='F_gsr.nii.gz',
                                    out_pf_name='pF_gsr.nii.gz',
                                    demean=True),
                            name='gsr_filter')

    # output identities
    rest_noise_output = setup_node(IdentityInterface(fields=out_fields,
                                                     mandatory_inputs=True),
                                    name="rest_noise_output")

    # Connect the nodes
    wf.connect([
                # tsnr
                (rest_noise_input, tsnr, [("in_file", "in_file")]),

                # artifact detection
                (rest_noise_input, art, [("in_file",       "realigned_files"),
                                         ("motion_params", "realignment_parameters"),
                                         ("brain_mask",    "mask_file"),
                                        ]),

                # calculte motion regressors
                (rest_noise_input,  motion_regs, [("motion_params", "motion_params")]),

                # create motion and confound regressors parameters file
                (art,           motart_pars, [("norm_files",    "comp_norm"),
                                              ("outlier_files", "outliers"),
                                             ]),
                (motion_regs,   motart_pars, [("out_files", "motion_params")]),

                # motion filtering
                (rest_noise_input, motion_filter, [("in_file", "in_file"),
                                                   (("in_file", rename, "_filtermotart"), "out_res_name"),
                                                  ]),
                (motart_pars,      motion_filter, [(("out_files", selectindex, 0), "design")]),

                # output
                (tsnr,             rest_noise_output, [("tsnr_file",          "tsnr_file")]),
                (motart_pars,      rest_noise_output, [("out_files",          "motion_regressors")]),
                (motion_filter,    rest_noise_output, [("out_res",            "motion_corrected")]),
                (art,              rest_noise_output, [("displacement_files", "art_displacement_files"),
                                                       ("intensity_files",    "art_intensity_files"),
                                                       ("norm_files",         "art_norm_files"),
                                                       ("outlier_files",      "art_outlier_files"),
                                                       ("plot_files",         "art_plot_files"),
                                                       ("statistic_files",    "art_statistic_files"),
                                                      ]),
                ])

    last_filter = motion_filter

    # compcor filter
    if filters['compcor_csf'] or filters['compcor_wm']:
        wf.connect([
                    # calculate compcor regressor and parameters file
                    (motart_pars,      compcor_pars,      [(("out_files", selectindex, 0),   "extra_regressors"),]),
                    (motion_filter,    compcor_pars,      [("out_res",                       "realigned_file"),]),

                    # the compcor filter
                    (motion_filter,    compcor_filter,    [("out_res",                        "in_file"),
                                                           (("out_res", rename, "_cleaned"),  "out_res_name"),
                                                          ]),
                    (compcor_pars,     compcor_filter,    [(("out_files", selectindex, 0),    "design")]),
                    #(compcor_pars,     compcor_filter,    [("components_file",  "design")]),
                    (rest_noise_input, compcor_filter,    [("brain_mask",                     "mask")]),

                    # output
                    (compcor_pars,     rest_noise_output, [("out_files",   "compcor_regressors")]),
                    #(compcor_pars,     rest_noise_output, [("components_file",   "compcor_regressors")]),
                    ])
        last_filter = compcor_filter

    # global signal regression
    if filters['gsr']:
        wf.connect([
            # calculate gsr regressors parameters file
            (last_filter,       gsr_pars, [("out_res",      "realigned_file")]),
            (rest_noise_input,  gsr_pars, [("brain_mask",   "mask_file")]),

            # the output file name
            (rest_noise_input, gsr_filter,  [("brain_mask", "mask")]),
            (last_filter,       gsr_filter, [("out_res", "in_file"),
                                             (("out_res", rename, "_gsr"), "out_res_name"),
                                            ]),
            (gsr_pars,          gsr_filter, [(("out_files", selectindex, 0), "design")]),

            # output
            (gsr_pars,   rest_noise_output, [("out_files", "gsr_regressors")]),
        ])
        last_filter = gsr_filter

    # connect the final nuisance correction output node
    wf.connect([(last_filter, rest_noise_output, [("out_res", "nuis_corrected")]),])

    if filters['compcor_csf'] and filters['compcor_wm']:
        mask_merge = setup_node(Merge(2), name="mask_merge")
        wf.connect([
                    ## the mask for the compcor filter
                    (rest_noise_input, mask_merge,   [("wm_mask",    "in1")]),
                    (rest_noise_input, mask_merge,   [("csf_mask",   "in2")]),
                    (mask_merge,       compcor_pars, [("out",        "mask_file")]),
                  ])

    elif filters['compcor_csf']:
        wf.connect([
                    ## the mask for the compcor filter
                    (rest_noise_input, compcor_pars, [("csf_mask",   "mask_file")]),
                  ])

    elif filters['compcor_wm']:
        wf.connect([
                    ## the mask for the compcor filter
                    (rest_noise_input, compcor_pars, [("wm_mask",   "mask_file")]),
                  ])

    return wf
