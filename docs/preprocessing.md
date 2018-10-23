
# Pre-processing pipelines
Here I detail each of the different pre-processing pipelines and their configuration parameters.

<!-- TOC depthFrom:2 depthTo:6 withLinks:1 updateOnSave:1 orderedList:0 -->

- [T1-weighted MRI](#t1-weighted-mri)
- [Resting-state fMRI (RS-fMRI)](#resting-state-fmri)
- [Diffusion MRI](#diffusion-mri)
- [Positron-Emission Tomography](#positron-emission-tomography)
- [T1-Weighted MRI & Positron Emission Tomography](#t1-weighted-mri-\&-positron-emission-tomography)

<!-- /TOC -->

## T1-weighted MRI
[<a href="https://github.com/Neurita/pypes/blob/master/docs/img/spm_anat_preproc_workflow.png?raw=true" target="_blank">graph</a>]
This pipeline will bias-field correct, segment the tissues, and register a
**T1-weighted image** to MNI.

It is based in **ANTS and SPM12**.
It is implemented in [`neuro_neuro_pypes.anat.attach_spm_anat_preprocessing`](https://github.com/Neurita/pypes/blob/master/neuro_pypes/anat/preproc.py).

A cortical thickness method is enabled with the `anat_preproc.do_cortical_thickness` boolean field.
This performs the **SPM+DiReCT** method described in [(Schwarz et al., 2016)](http://dx.doi.org/10.1016/j.nicl.2016.05.017),
which will use ANTs' KellyKapowski tool on SPM12 tissue segmentations.

These are the **steps**:

1. **Bias-field** correction using ANTS/N4BiasFieldCorrection.
2. Brain **tissue segmentation** and **spatial normalization** with SPM12 New Segment.
3. **Tissue maps spatial normalization** to MNI using SPM12 Normalize.
4. Create a **brain mask** based on tissue segmentations.

[optional]

5. **Warp atlas** (or any file in SPM12-MNI space) to anatomical space.
6. Measure **Cortical Thickness** with the SPM+DiReCT method.

##### Related settings

```yaml
spm_dir: "~/Software/matlab_tools/spm12"

# anatomical image pre-processing
normalize_atlas: True
atlas_file: '/home/hansel/data/std_brains/atlases/hammers/Hammers_mith_atlas_n30r83_SPM5.nii.gz'

# this is similar to the SPM+DiReCT method described here:
# http://dx.doi.org/10.1016/j.nicl.2016.05.017
anat_preproc.do_cortical_thickness: True

# these are the KellyKapowski default parameters
# also used in antsCorticalThickness
direct.convergence: "[45,0.0,10]"
direct.gradient_step: 0.025
direct.smoothing_variance: 1.0
direct.smoothing_velocity_field: 1.5
direct.use_bspline_smoothing: False
direct.number_integration_points: 10
direct.thickness_prior_estimate: 10
```


## Resting-state fMRI
[<a href="https://github.com/Neurita/pypes/blob/master/docs/img/spm_rest_preproc_workflow.png?raw=true" target="_blank">graph</a>]
This pipeline preprocess fMRI data for **resting-state fMRI** (rs-fMRI) analyses.
It depends on the T1-weighted preprocessing pipeline.

It is based on **SPM12, nipype ArtifactDetect and TSNR, Nipy motion correction, and nilearn.**

It consists on two parts:
1. fMRI data cleaning ( [`neuro_pypes.fmri.clean.attach_fmri_cleanup_wf`](https://github.com/Neurita/pypes/blob/master/neuro_pypes/fmri/clean.py)) and
2. warping and smoothing ([`neuro_pypes.fmri.warp.attach_spm_warp_fmri_wf`](https://github.com/Neurita/pypes/blob/master/neuro_pypes/fmri/warp.py)).

The connection of both parts is in
[`neuro_pypes.fmri.resting._attach_rest_preprocessing`](https://github.com/Neurita/pypes/blob/master/neuro_pypes/fmri/resting.py).

It's also possible to create a **group template** if you set that in the
configuration file.

**Steps:**

1. **Trim** the first 6 seconds from the fMRI data.
2. **Slice-time correction** based on SPM12 SliceTiming.
3. **Motion correction** with nipy.SpaceTimeRealigner.
4. **Tissue co-registration** from anatomical space to fMRI space.
5. **Nuisance corrections** including:
    1. time-course SNR (TSNR) estimation,
    1. artifact detection (nipype.rapidART),
    1. motion estimation and filtering,
    1. signal component regression from different tissues (nipype ACompCor), and
    1. global signal regression (GSR).
6. **Bandpass time filter**. Settings: `rest_input.lowpass_freq: 0.1`, and
`rest_input.highpass_freq: 0.01`.
7. Spatial **smoothing**. `smooth_fmri.fwhm: 8`

The **slice-time correction** requires information in the headers of the NifTI files about acquisition
slice-timing. NifTI files generated from most DICOM formats with a recent
version of `dcm2niix` should have the necessary information.

The **nuisance correction** steps can be enabled/disabled and configurable through
the configuration file.
**Note:** There is one thing that can't be easily modified: you need to
perform component regression for at least one tissue, e.g., CSF.

In the same way as for the MRI + FDG-PET pipeline (described below),
there is **2 ways for registration**. This is configured through the `registration.fmri2mni`
option.

#### If registration.fmri2mni: True
8. Cleaned-up versions of **fMRI are directly warped** to MNI using SPM12
Normalize.
9. Smooth these warped images, in the same ways as the non-warped data,
according to `smooth_fmri.fwhm: 8`.

#### If registration.fmri2mni: False
8. **Co-register fMRI** to anatomical space.
9. Apply the anat-to-MNI warp field to warp the cleaned-up versions of
the fMRI data to MNI.
10. Smooth these warped images, in the same ways as the non-warped data,
according to `smooth_fmri.fwhm: 8`.


##### Related settings

```yaml
registration.fmri2mni: True

# degree of b-spline used for rs-fmri interpolation
coreg_rest.write_interp: 3

## the last volume index to discard from the timeseries. default: 0
trim.begin_index: 5

# REST (COBRE DB)
# http://fcon_1000.projects.nitrc.org/indi/retro/cobre.html
# Rest scan:
# - collected in the Axial plane,
# - series ascending,
# - multi slice mode and
# - interleaved.
stc_input.slice_mode: alt_inc
stc_input.time_repetition: 2
#stc_input.num_slices: 33

# fMRI PREPROCESSING
# for any fmri warping, except group template creation (look below)
fmri_warp.write_voxel_sizes: [2, 2, 2]

# mid-process registration for group template creation
fmri_grptemplate_warp.write_voxel_sizes: [2, 2, 2]

# GROUP fMRI TEMPLATE
spm_epi_grouptemplate_smooth.fwhm: 8

# path to a common EPI template, if you don't want the average across subjects
spm_epi_grouptemplate_smooth.template_file: ""

# bandpass filter frequencies in Hz.
rest_input.lowpass_freq: 0.1 # the numerical upper bound
rest_input.highpass_freq: 0.01 # the numerical lower bound

# fwhm of smoothing kernel [mm]
smooth_fmri.fwhm: 8

## CompCor rsfMRI filters (at least compcor_csf should be True).
rest_filter.compcor_csf: True
rest_filter.compcor_wm: False
rest_filter.gsr: False

# filters parameters
## the corresponding filter must be enabled for these.

# motion regressors upto given order and derivative
# motion + d(motion)/dt + d2(motion)/dt2 (linear + quadratic)
motion_regressors.order: 0
motion_regressors.derivatives: 1

# number of polynomials to add to detrend
motart_parameters.detrend_poly: 2

# Compute TSNR on realigned data regressing polynomials up to order 2
tsnr.regress_poly: 2

# Threshold to use to detect motion-related outliers when composite motion is being used
detect_artifacts.use_differences: [True, False]
detect_artifacts.parameter_source: NiPy
detect_artifacts.mask_type: file
detect_artifacts.use_norm: True
detect_artifacts.zintensity_threshold: 3
detect_artifacts.norm_threshold: 1

# Number of principal components to calculate when running CompCor. 5 or 6 is recommended.
compcor_pars.num_components: 6

# Number of principal components to calculate when running Global Signal Regression. 1 is recommended.
gsr_pars.num_components: 1
```


## Diffusion MRI
[<a href="https://github.com/Neurita/pypes/blob/master/docs/img/spm_pet_anat_dti_camino_workflow.png?raw=true" target="_blank">graph</a>]
This pipeline performs **Diffusion MRI** (DTI) correction and pre-processing, tensor-fitting and tractography.

It is based on **FSL Eddy, dipy, and UCL Camino**.

**Steps:**

1. **Eddy-currents and motion correction** through FSL Eddy.
2. **Non-Local Means** from dipy for image de-noising with a Rician filter.
3. **Co-register the anatomical** image to diffusion space.
4. **Rotate the b-vecs** based on motion estimation.

[optional]

5. **Warp an atlas** to diffusion space.

The **Eddy-currents correction** needs certain fields in the NifTI file
to be able to create input parameters for Eddy.
Any file converted from modern DICOM to NifTI with a recent version of `dcm2nii` or `dcm2niix`
should work.

The atlas warping is needed if you'll perform **tractography**.


##### Related settings
```yaml
normalize_atlas: True
atlas_file: ''

# degree of b-spline used for interpolation
coreg_b0.write_interp: 3
nlmeans_denoise.N: 12 # number of channels in the head coil
```


## Positron-Emission Tomography
This is a **spatial normalization** pipeline for Positron-Emission Tomography (**PET**) images.
This workflow has showed good registration results on FDG and FDOPA PET images.

It is based on **SPM12**. You can find its source code in
[`neuro_pypes.pet.warp.attach_spm_pet_preprocessing`](https://github.com/Neurita/pypes/blob/master/neuro_pypes/pet/warp.py).

**Steps:**

1. **Spatially normalize** FDG-PET to MNI using SPM12 Normalize.

There is a **group template** option for PET: first a group template
is created, then all subjects are normalized to this
group template.

##### Related settings
```yaml
# GROUP PET TEMPLATE
spm_pet_grouptemplate_smooth.fwhm: 8
# path to a common PET template, if you don't want the average across subjects
spm_pet_grouptemplate.template_file: ""
```

## T1-Weighted MRI & Positron Emission Tomography
[<a href="https://github.com/Neurita/pypes/blob/master/docs/img/spm_anat_pet_preproc_workflow.png?raw=true" target="_blank">graph</a>]
This is a **partial volume correction (PVC) and spatial normalization pipeline for PET images**.

It is based on **PETPVC, nilearn and SPM12**.
It is implemented in [`neuro_pypes.pet.mrpet.attach_spm_mrpet_preprocessing`](https://github.com/Neurita/pypes/blob/master/neuro_pypes/pet/mrpet.py).

This pipeline **depends** on the anatomical preprocessing pipeline.
There is **2 different ways of co-registration**, you can configure that by
setting the `registration.anat2pet` boolean option to `True` or `False`.

#### If registration.anat2pet: True
1. **Co-register** anatomical and tissues to PET space.
2. **Partial volume effect correction** (PVC) with PETPVC in PET space.
3. **Directly normalize PET to MNI** with SPM12 Normalize.

#### If registration.anat2pet: False
1. **Co-register** PET to anatomical space.
2. **PVC** with PETPVC in anatomical space.
3. **Normalize PET to MNI** with SPM12 Normalize applying the
anatomical-to-MNI warp field.

[optional]

5. **Warp atlas** from anatomical to PET space.


The PVC is done based on tissue segmentations from the anatomical pipeline.


##### Related settings
```yaml
normalize_atlas: True
atlas_file: ''

registration.anat2pet: False

# GROUP PET TEMPLATE with MR co-registration
spm_mrpet_grouptemplate_smooth.fwhm: 8
spm_mrpet_grouptemplate.do_petpvc: True
# path to a common PET template, if you don't want the average across subjects
spm_mrpet_grouptemplate.template_file: ""

# PET PVC
pvc.pvc: RBV
pvc.fwhm_x: 4.3
pvc.fwhm_y: 4.3
pvc.fwhm_z: 4.3
```
