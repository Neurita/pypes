# Command Line

Pypes includes a CLI to help you with result check-up and analysis.

It is called: `nitap`.

To check what it is capable of run: `nitap -h`.

To indicate sets of files in the output folder, it uses `hansel.Crumb`
notation, see it in the examples below.

## `nitap plot`

Produces a report with plots of the resulting brain images.

You can plot the images in pairs of background and foreground, and 
it also allows you to choose the number of slices you want.

Example:

```
nitap plot --bg "/data/hansel/cobre/{sid}/{session}/anat.nii.gz" --fg "/data/hansel/cobre/{sid}/{session}/mni_in_anat_space.nii.gz"
```


## `nitap motion`

Produces a Excel spreadsheet with the result from RAPIDArt for all
pre-processed fMRI subjects.

The value of the crumb arguments in the input path will be also put 
in the spreadsheet.

Example:

```
nitap motion -i "/home/alexandre/data/nuk/out/{group}/{sid}/session_0/rest/artifact_stats" -o motion.xls
```