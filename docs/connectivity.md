
# Resting-state fMRI connectivity analysis
**This pipeline is under development.**

This pipeline warps an atlas file to the fMRI space, and then
perform the connectivity measures with the pre-processed rs-fMRI data.

This pipeline depends on the anatomical and rs-fMRI pre-processing pipelines.

There is already an interface almost done in [`neuro_pypes.interfaces.nilearn.connectivity`](https://github.com/Neurita/pypes/blob/master/neuro_pypes/interfaces/nilearn/connectivity.py) for this.

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
