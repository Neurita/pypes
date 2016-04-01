# -*- coding: utf-8 -*-
"""
Resting-state fMRI specific nuisance correction filtering workflow.
"""

import nipype.pipeline.engine       as pe
from   nipype.algorithms.rapidart   import ArtifactDetect
from   nipype.interfaces.utility    import Function, IdentityInterface, Merge
from   nipype.interfaces            import fsl

from   ..utils import (setup_node,
                       selectindex,
                       _get_params_for,
                       rename, )


def motion_regressors(motion_params, order=0, derivatives=1):
    """Compute motion regressors upto given order and derivative

    motion + d(motion)/dt + d2(motion)/dt2 (linear + quadratic)
    """
    import os
    import numpy as np
    from nipype.utils.filemanip import filename_to_list

    out_files = []
    for idx, filename in enumerate(filename_to_list(motion_params)):
        params = np.genfromtxt(filename)
        out_params = params
        for d in range(1, derivatives + 1):
            cparams = np.vstack((np.repeat(params[0, :][None, :], d, axis=0),
                                 params))
            out_params = np.hstack((out_params, np.diff(cparams, d, axis=0)))
        out_params2 = out_params
        for i in range(2, order + 1):
            out_params2 = np.hstack((out_params2, np.power(out_params, i)))
        filename = os.path.join(os.getcwd(), "motion_regressor%02d.txt" % idx)
        np.savetxt(filename, out_params2, fmt="%.10f")
        out_files.append(filename)
    return out_files


def create_regressors(motion_params, comp_norm, outliers, detrend_poly=None):
    """Builds a regressor set comprising motion parameters, composite norm and
    outliers.
    The outliers are added as a single time point column for each outlier

    Parameters
    ----------
    motion_params: a text file containing motion parameters and its derivatives
    comp_norm: a text file containing the composite norm
    outliers: a text file containing 0-based outlier indices
    detrend_poly: number of polynomials to add to detrend

    Returns
    -------
    components_file: a text file containing all the regressors
    """
    import os
    import numpy as np
    from nipype.utils.filemanip import filename_to_list
    from scipy.special import legendre

    out_files = []
    for idx, filename in enumerate(filename_to_list(motion_params)):
        params = np.genfromtxt(filename)
        norm_val = np.genfromtxt(filename_to_list(comp_norm)[idx])
        out_params = np.hstack((params, norm_val[:, None]))
        try:
            outlier_val = np.genfromtxt(filename_to_list(outliers)[idx])
        except IOError:
            outlier_val = np.empty((0))
        for index in np.atleast_1d(outlier_val):
            outlier_vector = np.zeros((out_params.shape[0], 1))
            outlier_vector[index] = 1
            out_params = np.hstack((out_params, outlier_vector))
        if detrend_poly:
            timepoints = out_params.shape[0]
            X = np.empty((timepoints, 0))
            for i in range(detrend_poly):
                X = np.hstack((X, legendre(
                    i + 1)(np.linspace(-1, 1, timepoints))[:, None]))
            out_params = np.hstack((out_params, X))
        filename = os.path.join(os.getcwd(), "filter_regressor_%02d.txt" % idx)
        np.savetxt(filename, out_params, fmt="%.10f")
        out_files.append(filename)
    return out_files


def extract_noise_components(realigned_file, mask_file, num_components=5,
                             extra_regressors=None):
    """Derive components most reflective of physiological noise
    Parameters
    ----------
    realigned_file: a 4D Nifti file containing realigned volumes
    mask_file: a 3D Nifti file containing white matter + ventricular masks
    num_components: number of components to use for noise decomposition
    extra_regressors: additional regressors to add
    Returns
    -------
    components_file: a text file containing the noise components
    """
    import os
    import nibabel as nb
    import numpy as np
    import scipy as sp
    from   nipype.utils.filemanip import filename_to_list

    imgseries = nb.load(realigned_file)
    components = None
    for filename in filename_to_list(mask_file):
        mask = nb.load(filename).get_data()
        if len(np.nonzero(mask > 0)[0]) == 0:
            continue
        voxel_timecourses = imgseries.get_data()[mask > 0]
        voxel_timecourses[np.isnan(np.sum(voxel_timecourses, axis=1)), :] = 0
        # remove mean and normalize by variance
        # voxel_timecourses.shape == [nvoxels, time]
        X = voxel_timecourses.T
        stdX = np.std(X, axis=0)
        stdX[stdX == 0] = 1.
        stdX[np.isnan(stdX)] = 1.
        stdX[np.isinf(stdX)] = 1.
        X = (X - np.mean(X, axis=0)) / stdX
        u, _, _ = sp.linalg.svd(X, full_matrices=False)
        if components is None:
            components = u[:, :num_components]
        else:
            components = np.hstack((components, u[:, :num_components]))
    if extra_regressors:
        regressors = np.genfromtxt(extra_regressors)
        components = np.hstack((components, regressors))
    components_file = os.path.join(os.getcwd(), 'noise_components.txt')

    np.savetxt(components_file, components, fmt="%.10f")
    return components_file


def rapidart_artifact_detection():
    art = ArtifactDetect()
    art.inputs.use_differences      = [True, True]
    art.inputs.use_norm             = True
    art.inputs.zintensity_threshold = 9
    art.inputs.mask_type            = 'spm_global'
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

    Nipype Outputs
    --------------
    rest_noise_output.motion_corrected
        The fMRI motion corrected image.

    rest_noise_output.nuis_corrected
        The resulting nuisance corrected image.
        This will be the same as 'motion_corrected' if compcor is disabled.

    rest_noise_output.motion_regressors
        Motion regressors file.

    rest_noise_output.compcor_regressors
        CompCor regressors file.

    Returns
    -------
    rm_nuisance_wf: nipype Workflow
    """

    # Create the workflow object
    wf = pe.Workflow(name=wf_name)

    # input identities
    rest_noise_input = setup_node(IdentityInterface(fields=["in_file",
                                                            "brain_mask",
                                                            "wm_mask",
                                                            "csf_mask",
                                                            "motion_params"],
                                                    mandatory_inputs=True),
                                  name="rest_noise_input")

    # get the settings for filters
    filters = _get_params_for('rest_filter')

    # Use :class:`nipype.algorithms.rapidart` to determine which of the
    # images in the functional series are outliers based on deviations in
    # intensity or movement.
    art = setup_node(rapidart_artifact_detection(), name="detect_artifact")

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

    # Noise confounds regressors
    compcor_pars = setup_node(Function(input_names=['realigned_file',
                                                    'mask_file',
                                                    'num_components',
                                                    'extra_regressors'],
                                         output_names=['out_files'],
                                         function=extract_noise_components,),
                              name='compcor_pars')

    compcor_filter = setup_node(fsl.GLM(out_f_name='F.nii.gz',
                                        out_pf_name='pF.nii.gz',
                                        demean=True),
                                name='compcor_filter')

    # output identities
    rest_noise_output = setup_node(IdentityInterface(fields=[
                                                             "compcor_corrected",
                                                             "motion_corrected",
                                                             "nuis_corrected",
                                                             "motion_regressors",
                                                             "compcor_regressors",
                                                            ],
                                                     mandatory_inputs=True),
                                    name="rest_noise_output")

    # Connect the nodes
    wf.connect([
                # artifact detection
                (rest_noise_input, art, [("in_file",        "realigned_files"),
                                         ("motion_params",  "realignment_parameters"),
                                        ]),

                # calculte motion regressors
                (rest_noise_input,  motion_regs, [("motion_params", "motion_params")]),

                # create motion and confound regressors parameters file
                (art,           motart_pars, [("norm_files",    "comp_norm"),
                                              ("outlier_files", "outliers"),
                                             ]),
                (motion_regs,   motart_pars, [("out_files",     "motion_params")]),

                # motion filtering
                (rest_noise_input, motion_filter,   [("in_file",                            "in_file"),
                                                     (("in_file", rename, "_filtermotart"), "out_res_name"),
                                                    ]),
                (motart_pars,      motion_filter,   [(("out_files", selectindex, [0]),      "design")]),

                # output
                (motart_pars,      rest_noise_output, [("out_files",   "motion_regressors")]),
                (motion_filter,    rest_noise_output, [("out_res",     "motion_corrected")]),
              ])

    if filters['compcor_csf'] or filters['compcor_wm']:
        wf.connect([
                    # calculate compcor regressor and parameters file
                    (motart_pars,      compcor_pars,      [(("out_files", selectindex, [0]), "extra_regressors"),]),
                    (motion_filter,    compcor_pars,      [("out_res",                       "realigned_file"),]),

                    # the compcor filter
                    (motion_filter,    compcor_filter,    [("out_res",                        "in_file"),
                                                           (("out_res", rename, "_cleaned"),  "out_res_name"),
                                                          ]),
                    (compcor_pars,     compcor_filter,    [(("out_files", selectindex, [0]),  "design")]),
                    (rest_noise_input, compcor_filter,    [("brain_mask",                     "mask")]),

                    # output
                    (compcor_pars,     rest_noise_output, [("out_files",   "compcor_regressors")]),
                    (compcor_filter,   rest_noise_output, [("out_res",     "compcor_corrected")]),
                    (compcor_filter,   rest_noise_output, [("out_res",     "nuis_corrected")]),
                  ])
    else:
        wf.connect([
                    ## the motion corrected is the nuis_corrected result
                    (motion_filter, rest_noise_output, [("out_res", "nuis_corrected")]),
                  ])

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
