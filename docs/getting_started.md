# Getting started

In this document I show how to build and run neuroimaging pipelines using [Nipype](http://nipype.readthedocs.io) and [pypes](https://github.com/Neurita/pypes).


<!-- TOC depthFrom:2 depthTo:6 withLinks:1 updateOnSave:1 orderedList:0 -->

- [Preparing the data](#preparing-the-data)
- [Declaring the data](#declaring-the-data)
- [Setting up the pipelines](#setting-up-the-pipelines)
- [Setting pipeline parameters](#setting-pipeline-parameters)
- [Building the pipeline](#building-the-pipeline)
- [Running the pipelines](#running-the-pipelines)

<!-- /TOC -->

## Preparing the data
I recommend organizing your dataset first. [This guide](http://miykael.github.io/nipype-beginner-s-guide/prepareData.html) explains good ways to organize it. Usually I recommend the following tree: `{base_dir}/raw/{subject_id}/{session_id}/{image_files}`.

## Declaring the data
Now let's start programming. The first step is to declare and explain the file tree for the pipeline input.

Some imports first:

```python
import os.path as path

from hansel import Crumb
```

Let's imagine the base folder where our data is:

```python
base_dir = "/home/pyper/data/cobre/raw"
```

Now we specify the tree data structure using `hansel.Crumb`.
We are going to use `path.join` to build the path until the input files of the pipeline. Each part (subfolder) in this path that has many values, you enclose it with curly braces and put a sensible name for it. For example, in the level where all the folder named as subject identifiers, I will put `"{subject_id}"`.

```python
data_path = path.join(base_dir, "{subject_id}", "session_1", "{modality}", "{image}")

data_crumb = Crumb(data_path, ignore_list=[".*"])
#
```

Now in this way we have the data structure tree defined.
```python
print(data_crumb)
>>> Crumb("/home/pyper/data/cobre/raw/{subject_id}/session_1/{modality}/{image}")
```

## Setting up the pipelines

Once the data is defined, we have to specify the pipelines we want to run against these data. Let's imagine we want to run the anatomical and the function pre-processing pipelines.

First we import the attach functions to build those pipelines:

```python
from pypes.anat import attach_spm_anat_preprocessing
from pypes.fmri import attach_rest_preprocessing
```

Then we specify and give names to the pipelines:

```python
attach_functions = {"spm_anat_preproc": attach_spm_anat_preprocessing,
                    "spm_rest_preproc": attach_rest_preprocessing,}

```

One last thing, we have to link each file type to the corresponding pipeline. The `attach_spm_anat_preprocessing` expects an `anat` input and the `attach_rest_preprocessing` expects a `rest` input. So we have to specify the parameters in the data tree for each of these (check the argument names we used before):

```python
crumb_arguments = {'anat': [('modality', 'anat_1'),
                            ('image',    'mprage.nii.gz')],
                   'rest': [('modality', 'rest_1'),
                            ('image',    'rest.nii.gz')],
                  }
```

## Setting pipeline parameters

You can either have a separate configuration file (in JSON, YAML, .ini, etc...), dealt by [Kaptan](https://github.com/emre/kaptan). And you can also set configuration parameters with a Python `dict`.
Pypes provides a function named `update_config` where you can input a `dict` or a file path to update the global configuration state before building the pipelines.
I recommend YAML, have a look at the [example config file](pypes_config.yml).

```python
from pypes.config import update_config

if config_file_path:
    update_config(config_file_path)

if params_dict:
    update_config(params_dict)
```

## Building the pipeline

Finally we define the output and working folders:

```python
output_dir = path.join(path.dirname(base_dir), "out")

cache_dir = path.join(path.dirname(base_dir), "wd")
```

And we instantiate the workflow with the variables we created before:

```python
from pypes.io import build_crumb_workflow

wf = build_crumb_workflow(attach_funcs,
                          data_crumb=data_crumb,
                          in_out_kwargs=crumb_arguments,
                          output_dir=output_dir,
                          cache_dir=cache_dir,)
```

## Running the pipelines

To run the workflow we can use two functions, with or without debugging.
The one in debug mode will call a debugger if any exception is raised in the end. These functions are: `pypes.run.run_debug` and `pypes.run.run_wf`.
They have the same arguments as the workflow `wf.run()` function.

```python
from pypes.run import run_debug

run_debug(wf, plugin="MultiProc", n_cpus=4)
```
