# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:

import operator

import nipype.pipeline.engine    as pe
from   nipype.algorithms.misc    import TSNR
from nipype.algorithms.rapidart  import ArtifactDetect
from   nipype.interfaces.utility import Function, IdentityInterface
from   nipype.interfaces         import fsl

from   ..utils import (setup_node,
                       _get_params_for,)



def bandpass_filter(files, lowpass_freq, highpass_freq, tr=2):
    """Bandpass filter the input files

    Parameters
    ----------
    files: list of 4d nifti files
    lowpass_freq: cutoff frequency for the low pass filter (in Hz)
    highpass_freq: cutoff frequency for the high pass filter (in Hz)
    tr: float
        The repetition time in seconds. The inverse of sampling rate (in Hz)
    """
    import os
    import nibabel as nb
    import numpy as np
    from nipype.utils.filemanip import (filename_to_list,
                                        list_to_filename,
                                        split_filename)

    fs = 1/tr

    out_files = []
    for filename in filename_to_list(files):
        path, name, ext = split_filename(filename)
        out_file = os.path.join(os.getcwd(), name + '_bp' + ext)
        img = nb.load(filename)
        timepoints = img.shape[-1]
        F = np.zeros((timepoints))
        lowidx = int(timepoints / 2) + 1
        if lowpass_freq > 0:
            lowidx = np.round(float(lowpass_freq) / fs * timepoints)
        highidx = 0
        if highpass_freq > 0:
            highidx = np.round(float(highpass_freq) / fs * timepoints)
        F[highidx:lowidx] = 1
        F = ((F + F[::-1]) > 0).astype(int)
        data = img.get_data()
        if np.all(F == 1):
            filtered_data = data
        else:
            filtered_data = np.real(np.fft.ifftn(np.fft.fftn(data) * F))
        img_out = nb.Nifti1Image(filtered_data, img.affine, img.header)
        img_out.to_filename(out_file)
        out_files.append(out_file)
    return list_to_filename(out_files)


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


def build_filter1(motion_params, comp_norm, outliers, detrend_poly=None):
    """Builds a regressor set comprisong motion parameters, composite norm and
    outliers

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
        filename = os.path.join(os.getcwd(), "filter_regressor%02d.txt" % idx)
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
    from nipype.utils.filemanip import filename_to_list

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

def trash():
    """Use :class:`nipype.algorithms.rapidart` to determine which of the
    images in the functional series are outliers based on deviations in
    intensity or movement.
    """
    art = setup_node(interface=ArtifactDetect(), name="detect_artifact")
    art.inputs.use_differences = [True, True]
    art.inputs.use_norm = True
    art.inputs.zintensity_threshold = 9
    art.inputs.mask_type = 'spm_global'
    art.inputs.parameter_source = 'NiPy'

    """Here we are connecting all the nodes together. Notice that we add the merge node only if you choose
    to use 4D. Also `get_vox_dims` function is passed along the input volume of normalise to set the optimal
    voxel sizes.
    """

    wf.connect([(name_unique, realign, [('out_file', 'in_file')]),
                (realign, art, [('out_file', 'realigned_files')]),
                (realign, art, [('par_file', 'realignment_parameters')]),
                ])

    def selectindex(files, idx):
        import numpy as np
        from nipype.utils.filemanip import filename_to_list, list_to_filename
        return list_to_filename(np.array(filename_to_list(files))[idx].tolist())

    mask = Node(fsl.BET(), name='getmask')
    mask.inputs.mask = True
    wf.connect(calc_median, 'median_file', mask, 'in_file')
    # get segmentation in normalized functional space

    def merge_files(in1, in2):
        out_files = filename_to_list(in1)
        out_files.extend(filename_to_list(in2))
        return out_files


def rest_noise_filter_wf(wf_name='rest_noise_removal'):
    """ Create a resting-state fMRI noise removal node.

    Nipype Inputs
    -------------


    Nipype Outputs
    --------------


    Returns
    -------
    rmnoise_node: nipype.Node

    """

    # Create the workflow object
    wf = pe.Workflow(name=wf_name)


    # input identities
    rest_noise_input = setup_node(IdentityInterface(fields=["in_file",
                                                            "wm_mask",
                                                            "csf_mask",
                                                            "motion_params"], mandatory_inputs=True),
                                  name="rest_noise_input")

    #filters = _get_params_for('rest_filter')
    #sorted_filters = sorted(filters.items(), key=operator.itemgetter(1))
    #import ipdb; ipdb.set_trace()

    # Use :class:`nipype.algorithms.rapidart` to determine which of the
    # images in the functional series are outliers based on deviations in
    # intensity or movement.
    art = setup_node(interface=ArtifactDetect(), name="detect_artifact")
    art.inputs.use_differences = [True, True]
    art.inputs.use_norm = True
    art.inputs.zintensity_threshold = 9
    art.inputs.mask_type = 'spm_global'
    art.inputs.parameter_source = 'NiPy'

    # Compute motion regressors
    motion_reg = setup_node(Function(input_names=['motion_params',
                                                  'order',
                                                  'derivatives'],
                                     output_names=['out_files'],
                                     function=motion_regressors,),
                            name='motion_regression')

    # Create a filter to remove motion and art confounds
    motion_filt = setup_node(Function(input_names=['motion_params',
                                                   'comp_norm',
                                                   'outliers',
                                                   'detrend_poly'],
                                      output_names=['out_files'],
                                      function=build_filter1),
                             name='motion_filter')


    wf.connect(art, 'norm_files', createfilter1, 'comp_norm')
    wf.connect(art, 'outlier_files', createfilter1, 'outliers')

    filter1 = setup_node(fsl.GLM(out_f_name='F_mcart.nii.gz',
                              out_pf_name='pF_mcart.nii.gz',
                              demean=True),
                      iterfield=['in_file', 'design', 'out_res_name'],
                      name='filtermotion')

    wf.connect(realign, 'out_file', filter1, 'in_file')
    wf.connect(realign, ('out_file', rename, '_filtermotart'),
               filter1, 'out_res_name')
    wf.connect(createfilter1, 'out_files', filter1, 'design')

    createfilter2 = setup_node(Function(input_names=['realigned_file', 'mask_file',
                                                  'num_components',
                                                  'extra_regressors'],
                                     output_names=['out_files'],
                                     function=extract_noise_components,
                                     imports=imports),
                            iterfield=['realigned_file', 'extra_regressors'],
                            name='makecompcorrfilter')
    createfilter2.inputs.num_components = num_components

    wf.connect(createfilter1, 'out_files', createfilter2, 'extra_regressors')
    wf.connect(filter1, 'out_res', createfilter2, 'realigned_file')
    wf.connect(registration, ('outputspec.segmentation_files', selectindex, [0, 2]),
               createfilter2, 'mask_file')

    filter2 = setup_node(fsl.GLM(out_f_name='F.nii.gz',
                              out_pf_name='pF.nii.gz',
                              demean=True),
                      iterfield=['in_file', 'design', 'out_res_name'],
                      name='filter_noise_nosmooth')
    wf.connect(filter1, 'out_res', filter2, 'in_file')
    wf.connect(filter1, ('out_res', rename, '_cleaned'),
               filter2, 'out_res_name')
    wf.connect(createfilter2, 'out_files', filter2, 'design')
    wf.connect(mask, 'mask_file', filter2, 'mask')






    # output identities
    rest_noise_output = setup_node(IdentityInterface(fields=[
                                                             "motion_corrected",
                                                             "motion_params",
                                                             "tissues",
                                                             "anat",
                                                             "time_filtered",
                                                             "brain_mask",
                                                            ],
                                                      mandatory_inputs=True),
                                    name="rest_noise_output")


    # Connect the nodes
    wf.connect([
                # trim
                (rest_noise_input, motion_reg,     [("motion_params",   "motion_params")]),
                (motion_reg,       motion_filt,    [("out_files",       "motion_params")]),

    wf.connect(motreg, 'out_files', createfilter1, 'motion_params')
                (bandpass_filter, rest_output, [("out_file",                  "time_filtered")]),
              ])

    return wf
