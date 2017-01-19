# Datasets

Pypes already has a few sets of workflows which you can use. You would only have
to download the raw dataset or organize your data in the same way explained below,
and setup the [pipeline parameters](pypes_config.yml).

<!-- TOC depthFrom:2 depthTo:6 withLinks:1 updateOnSave:1 orderedList:0 -->

- [COBRE](#COBRE)
- [A clinical dataset](#a-clinical-dataset)

<!-- /TOC -->


## COBRE
The [COBRE dataset](http://fcon_1000.projects.nitrc.org/indi/retro/cobre.html)
consists of raw anatomical and functional MR data from 72 patients with Schizophrenia
and 75 healthy controls.

Once you [download](http://fcon_1000.projects.nitrc.org/indi/retro/cobre.html)
this dataset, you will find a file structure tree as this:
`{base_dir}/cobre/{subject_id}/session_1/{modality}/{image}`.

In pypes we provide a function readily prepared to run anatomical and
resting-state fMRI preprocessing on this database:
[`pypes.datasets.cobre_crumb_workflow`](https://github.com/Neurita/pypes/blob/master/pypes/datasets.py).

### How to use it

```python
import os.path as path

from hansel import Crumb
from pypes.datasets import cobre_crumb_workflow
from pypes.run import run_debug

# we downloaded the database in:
base_dir = '/home/pyper/data/cobre/raw'
cobre_tree = path.join('{subject_id}', 'session_1', '{modality}', '{image}')

# we define the database tree
cobre_crumb = Crumb(path.join(base_dir, cobre_tree), ignore_list=['.*'])

# output and working dir
output_dir = path.join(path.dirname(base_dir), 'out')
cache_dir  = path.join(path.dirname(base_dir), 'wd')

# we have a configuration file in:
config_file = path.join(path.dirname(base_dir), 'pypes_config.yml')

# we choose what pipeline set we want to run.
# the choices are: 'spm_anat_preproc', 'spm_rest_preproc'
wf_name = 'spm_rest_preproc' # for MPRAGE and rs-fMRI preprocessing

# instantiate the workflow
wf = cobre_crumb_workflow(wf_name     = wf_name,
                          data_crumb  = cobre_crumb,
                          cache_dir   = cache_dir,
                          output_dir  = output_dir,
                          config_file = config_file,
                         )

# run it
run_debug(wf, plugin='MultiProc', n_cpus=4)
```


## My clinical dataset

Sadly, this is not public available.

The dataset we are working in our department has a very similar folder structure
as COBRE:  `{base_dir}/{subject_id}/{session_id}/{image}`.
If you organize your data in the same way, you can directly use the function
[`pypes.datasets.clinical_crumb_workflow`](https://github.com/Neurita/pypes/blob/master/pypes/datasets.py).
Have a look at the [_clinical_wf_setup](https://github.com/Neurita/pypes/blob/master/pypes/datasets.py)
function to see all the sets of pipelines you can pick.


```python
import os.path as path

from hansel import Crumb
from pypes.datasets import clinical_crumb_workflow
from pypes.run import run_debug

# we downloaded the database in:
base_dir = '/home/pyper/data/nuk/raw'
data_tree = path.join('{subject_id}', '{session_id}', '{image}')

# we define the database tree
data_crumb = Crumb(path.join(base_dir, data_tree), ignore_list=['.*'])

# output and working dir
output_dir = path.join(path.dirname(base_dir), 'out')
cache_dir  = path.join(path.dirname(base_dir), 'wd')

# we have a configuration file in:
config_file = path.join(path.dirname(base_dir), 'pypes_config.yml')

# we choose what pipeline set we want to run.
# there are many choices
wf_name = 'spm_anat_pet_tpm_pvc' # MPRAGE preprocessing, PET MNI group template, PET PVC, and PET normalization to group template
# another could be 'anat_dti_camino' for MPRAGE and DTI/tractography

# instantiate the workflow
wf = clinical_crumb_workflow(wf_name     = wf_name,
                             data_crumb  = data_crumb,
                             cache_dir   = cache_dir,
                             output_dir  = output_dir,
                             config_file = config_file,
                             )

# run it
run_debug(wf, plugin='MultiProc', n_cpus=4)
```
