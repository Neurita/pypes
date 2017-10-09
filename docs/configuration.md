
# Configuration

The main nodes in the pipelines are configurable through a
configuration file.

We recommend using the YAML (https://en.wikipedia.org/wiki/YAML)
format for this file.

To change the default value of a node parameter you
have to add to the configuration file an entry for the value
you want, as in: `<node_name>.<parameter_name>: <value>`

For example, let's say we have a
[`spm.Normalize12`](http://nipype.readthedocs.io/en/latest/interfaces/generated/nipype.interfaces.spm.preprocess.html#normalize12)
node named `anat_warp` in one of the workflows.
We want to set the value of the parameter `bias_regularization` to `0.1`.
The entry should be:

```
anat_warp.bias_regularization: 0.1
```

The most direct way of checking the names of the nodes defined in
`Pypes`, is by looking into the source code.
You can also have a look at the
[default configuration file](https://github.com/Neurita/pypes/blob/master/docs/pypes_config.yml)
where we list the settings for the main nodes in each workflow.
In addition we have added a **"Related settings"** section for each
workflow description where we list the most important configuration
settings, taken from this file.
