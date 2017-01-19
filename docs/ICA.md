
# fMRI Independent Component Analysis (ICA)
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

### See also
- [Plotting ICA results](plotting.md#ICA)
