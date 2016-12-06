
# Description of the pipelines
Here I explain the key points of each of the post-processing pipelines.
I also list the configuration parameters relevant to each of them.


## fMRI Independent Component Analysis (ICA)
This pipeline performs ICA on fMRI images. It is based on nilearn, and you
can choose between CanICA and DictLearning.
There is one version for one functional image, in `attach_canica` and
another
for a group ICA (GICA) in `attach_concat_canica`. However, probably the GICA
approach should be further tested on real data.

It depends on the RS-fMRI pipeline.
This is implemented in
[`pypes.postproc.decompose`](https://github.com/Neurita/pypes/blob/master/pypes/postproc/decompose.py).


##### Related settings
```yaml
# INDEPENDENT COMPONENTS ANALYSIS
## True to perform CanICA
rest_preproc.canica: False

# CanICA settings
canica.algorithm: 'canica' # choices: 'canica', 'dictlearning'
canica.do_cca: True
canica.standardize: True
canica.n_components: 20
canica.threshold: 2.0
canica.smoothing_fwhm: 8
#canica.random_state: 0

canica_extra.plot: True
canica_extra.plot_thr: 2.0 # used if threshold is not set in the ICA
```

## RS-fMRI Connectivity
**This pipeline is under development.**

This pipeline would need to warp an atlas file with the fMRI image, and then
perform the connectivity measures.
There is already an interface almost done in [`pypes.interfaces.nilearn.connectivity`](https://github.com/Neurita/pypes/blob/master/pypes/interfaces/nilearn/connectivity.py) for this.

##### Related settings
```yaml
normalize_atlas: True
atlas_file: ''

# RS-fMRI CONNECTIVITY
## if atlas_file is defined, perform connectivity analysis
rest_preproc.connectivity: True
## if further smoothing (remember the output of the rest workflow is already smoothed)
rest_connectivity.standardize: False
rest_connectivity.kind: correlation # choices: "correlation", "partial correlation", "tangent", "covariance", "precision".
rest_connectivity.smoothing_fwhm: 8
#rest_connectivity.resampling_target: # choices: "mask", "maps" or undefined.
rest_connectivity.atlas_type: labels # choices: "labels", "probabilistic".
```


## Diffusion Tractography
[<a href="https://github.com/Neurita/pypes/blob/master/docs/img/spm_pet_anat_dti_camino_workflow.png?raw=true" target="_blank">graph</a>]
This pipeline performs DTI tensor model fitting and tractography.

It is based on UCL Camino.
It depends on the MPRAGE and the DTI pipeline.

1. DTI fit.
2. ROIxROI atlas-based deterministic tractography.
3. Connectivity matrices: one with the number of tracts for each pair of
ROIs, the other with average tract FA values for each pair.

##### Related settings
```yaml
normalize_atlas: True
atlas_file: ''

# Camino Tractography
track.curvethresh: 50
track.anisthresh: 0.2
```
