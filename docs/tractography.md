
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
