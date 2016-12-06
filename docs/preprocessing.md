
# Description of the pipelines
Here I explain the key points of each of the different preprocessing pipelines.
I also list the configuration parameters relevant to each of them.

## Anatomical MRI (MPRAGE)
[<a href="https://github.com/Neurita/pypes/blob/master/docs/img/spm_anat_preproc_workflow.png?raw=true" target="_blank">graph</a>]
This pipeline will bias-field correct, segment the tissues, and register a T1-weighted image to MNI.

It is based in ANTS and SPM12.
It is implemented in [`pypes.anat.attach_spm_anat_preprocessing`](https://github.com/Neurita/pypes/blob/master/pypes/anat.py).

These are the steps:

1. Bias-field correction using ANTS/N4BiasFieldCorrection.
2. Brain tissue segmentation and spatial normalization with SPM12 New Segment.
3. Spatial normalization of tissue maps to MNI using SPM12 Normalize.
4. Create a brain mask based on tissue segmentations.

[optional]

5. Warp atlas (or any file in SPM12-MNI space) to anatomical space.

##### Related settings

```yaml
spm_dir: "~/Software/matlab_tools/spm12"

normalize_atlas: True
atlas_file: ''
```

## FDG-PET
This is a spatial normalization pipeline for FDG-PET images. Here I say specifically FDG-PET because I haven't tested
this too much for other PET tracers.

It is based on SPM12.
It is implemented in [`pypes.pet.warp.attach_spm_pet_preprocessing`](https://github.com/Neurita/pypes/blob/master/pypes/pet/warp.py).

1. Use SPM12 Normalize to spatially normalize FDG-PET to MNI.

There is a group-template option of this: first a group template
is created, then all FDG-PET are images are normalized to this
group template.

##### Related settings
```yaml
# GROUP PET TEMPLATE
spm_pet_grouptemplate_smooth.fwhm: 8
# path to a common PET template, if you don't want the average across subjects
spm_pet_grouptemplate.template_file: ""
```

## MPRAGE + FDG-PET
[<a href="https://github.com/Neurita/pypes/blob/master/docs/img/spm_anat_pet_preproc_workflow.png?raw=true" target="_blank">graph</a>]
This is a partial volume correction and spatial normalization pipeline
for FDG-PET images.

It is based on PETPVC, nilearn and SPM12.
It is implemented in [`pypes.pet.mrpet.attach_spm_mrpet_preprocessing`](https://github.com/Neurita/pypes/blob/master/pypes/pet/mrpet.py).

This pipeline depends on the anatomical preprocessing pipeline.
There is 2 ways of doing the co-registration, you can configure that by
setting the `registration.anat2pet` boolean option to `True` or `False`.

#### If registration.anat2pet: True
1. Co-register anatomical and tissues to PET space.
2. Partial volume effect correction (PVC) with PETPVC in PET space.
This is done based on tissue segmentations from the anatomical pipeline.
3. Use SPM12 Normalize to normalize FDG-PET to MNI.

#### If registration.anat2pet: False
1. Co-register FDG-PET to anatomical space.
2. PVC with PETPVC in anatomical space.
3. Normalize PET to MNI with SPM12 Normalize applying the
anatomical-to-MNI warp field.

[optional]

5. Warp atlas from anatomical to PET space.

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
rbvpvc.pvc: RBV
rbvpvc.fwhm_x: 4.3
rbvpvc.fwhm_y: 4.3
rbvpvc.fwhm_z: 4.3
```

## Resting-state fMRI (RS-fMRI)
[<a href="https://github.com/Neurita/pypes/blob/master/docs/img/spm_rest_preproc_workflow.png?raw=true" target="_blank">graph</a>]
This pipeline preprocess fMRI data for resting-state fMRI analyses.
It depends on the MPRAGE preprocessing pipeline.

It is based on SPM12, nipype ArtifactDetect and TSNR, Nipy motion correction,
and nilearn.
It consists on two parts, the first is for data cleaning and the second for
warping and smoothing. The first is implemented in
[`pypes.fmri.clean.attach_fmri_cleanup_wf`](https://github.com/Neurita/pypes/blob/master/pypes/fmri/clean.py) and the latter is implemented in
[`pypes.fmri.warp.attach_spm_warp_fmri_wf`](https://github.com/Neurita/pypes/blob/master/pypes/fmri/warp.py). Both parts are connected and
are also prepared for create a common group template if you set that in the
configuration file. This connection is implemented in
[`pypes.fmri.rest._attach_rest_preprocessing`](https://github.com/Neurita/pypes/blob/master/pypes/fmri/rest.py).

1. Trim the first 6 seconds from the fMRI data.
2. Slice-time correction based on SPM12 SliceTiming.
This requires information in the headers of the files about acquisition
slice-timing. NifTI files generated from most DICOM formats with a recent
version of `dcm2niix` should have the necessary information.
3. Motion correction with nipy.SpaceTimeRealigner.
4. Co-registration of tissues in anatomical space to fMRI space.
5. Nuisance correction including time-course SNR (TSNR) estimation,
artifact detection (nipype.rapidART), motion estimation and filtering, signal
component regression from different tissues (nipype ACompCor) and global
signal regression (GSR).
These are configurable through the configuration file.
There is one thing that can't be easily modified is that, you need to
perform component regression for at least one tissue, e.g., CSF.
6. Bandpass time filter. Settings: `rest_input.lowpass_freq: 0.1`, and
`rest_input.highpass_freq: 0.01`.
7. Spatial smoothing. `smooth_fmri.fwhm: 8`

In the same way as for the MRI + FDG-PET pipeline, there is 2 ways for
registration. This is configured through the `registration.anat2fmri`
option.

#### If registration.anat2fmri: True
8. Cleaned-up versions of fMRI are directly warped to MNI using SPM12
Normalize.
9. Smooth these warped images, in the same ways as the non-warped data,
according to `smooth_fmri.fwhm: 8`.

#### If registration.anat2fmri: False
8. Co-register fMRI to anatomical space.
9. Apply the anat-to-MNI warp field to warp the cleaned-up versions of
the fMRI data to MNI.
10. Smooth these warped images, in the same ways as the non-warped data,
according to `smooth_fmri.fwhm: 8`.


##### Related settings

```yaml
registration.anat2fmri: True

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
fmri_warp.write_voxel_sizes: [2, 2, 2]
fmri_grptemplate_warp.write_voxel_sizes: [2, 2, 2]

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


## Diffusion MRI (DTI)
[<a href="https://github.com/Neurita/pypes/blob/master/docs/img/spm_pet_anat_dti_camino_workflow.png?raw=true" target="_blank">graph</a>]
This pipeline performs Diffusion MRI correction and pre-processing., tensor-fitting and tractography
it is based on FSL Eddy, dipy, and UCL Camino.

1. Eddy currents and motion correction through FSL Eddy.
It needs certain fields in the NifTI file to be able to create input
parameters for Eddy. Any file converted from modern DICOM to NifTI with a
recent version of `dcm2nii` or `dcm2niix` should work.
2. Non-Local Means from dipy for image de-noising with a Rician filter.
3. Co-register the anatomical image to diffusion space.
4. Rotate the b-vecs based on motion estimation.

[optional]

5. Warp an atlas to diffusion space (needed if you'll perform tractography).

##### Related settings
```yaml
normalize_atlas: True
atlas_file: ''

# degree of b-spline used for interpolation
coreg_b0.write_interp: 3
nlmeans_denoise.N: 12 # number of channels in the head coil
```
