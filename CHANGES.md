Changelog
=========

Version 1.1.2 (16.04.2018)
--------------------------

### Features:

- Add `--extra` parameter to `nitap motion` for extra subject metadata.
- Bump nipype version to 1.0.2. Few fixes included.


Version 1.1.1 (08.03.2018)
--------------------------

### Bugfixes:

- Fix wrong import in cli/cli.py file.
- Fix bugs in `nitap plot` command.


Version 1.1.0 (01.03.2018)
--------------------------

### Features:

- Add `nitap` CLI with `plot` and `motion` subcommands.


### Docs:

- Added a new section for the command line interface.


### Code changes:

- Changed all import paths to absolute paths.
- Line width changed to 120.
- Many code style changes towards automatised lint checks.


Version 1.0.1
-------------

- Bug fix: some wrong import names with the module name change.
- Bug fix: added missing requirements.


Version 1.0.0
-------------

- API break: Change module name from `pypes` to `neuro_pypes`.
- Documentation update.


Version 0.3.5
-------------

- Fix connection naming in the fMRI pipeline atlas output.


Version 0.3.4
-------------

- Fix connection naming in the fMRI pipeline atlas output.


Version 0.3.3
-------------

- Fix output renaming in SPM+DiReCT pipeline.

- Add `black_bg` parameter to multi slices plotting functions.

- Add motion_stats_sheet to fmri utils.

- General fix on import ordering which was raising errors from new Nipype.

Version 0.3.2
-------------

- Add atlas warping to fMRI pipeline.

- Fixes in fMRI pipeline.

Version 0.3 - 0.3.1
-------------------
- Add motion statistics measures in DTI pipeline.

- Improve documentation.

- Add SPM+DiReCT option to anatomical preprocessing.

- Add `invoke` examples in tasks.py.

- Style and bug fixes.

- Fix imports for latest `nipype` version.

- Add `plot_ortho_slices` function to `nilearn` interface.


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
