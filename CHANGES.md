Changelog
=========

Version 0.3 - 0.3.2
-------------------
- Add motion statistics measures in DTI pipeline.

- Improve documentation.

- Add ANTs' KellyKapowski interface.

- Add SPM+DiReCT option to anatomical preprocessing.

- Add `invoke` examples in tasks.py.

- Style and bug fixes.

- Fix imports for latest `nipype` version.

- Add `plot_ortho_slices` function to `nilearn` interface.

- Add anatomical co-registered files to the DTI pre-processing
workflow output.

- Add `motion_stats_sheet` function to `fmri.utils`.

- Refactor PET/MR pre-processing.

- Add options: `mrpet.do_tissue_pvc` and `mrpet.do_gm_normalization`.

- Set `fmri_grptemplate_warp.write_voxel_sizes` and `fmri_warp.write_voxel_sizes`
configuration settings to [2, 2, 2].


Version 0.2.0 - 0.2.1
---------------------
- Add mkdocs documentation and readthedocs setup.

- Add registration options for PET and fMRI. Add options for this in config.

- Add group template for fMRI and PET pipelines.

- Add brain_mask output to anatomical MR pipeline.

- Rename/re-structure fMRI pipeline and files.

- Improve PET output file namings.

- Add settings for Camino in config.

- Fix bug in _get_n_slices in slicetime_params.

- Bug fixes.
