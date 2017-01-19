
# Configuration

The main nodes in the pipelines are configurable through a configuration file. 
To change a default value of configuration parameter of a node you have to add 
to the configuration file an entry for the value you want. 

For example, let's say we have a [`spm.Normalize12`](http://nipype.readthedocs.io/en/latest/interfaces/generated/nipype.interfaces.spm.preprocess.html#normalize12) 
node named `anat_warp` in one of the workflows. We want to set the value of the parameter
`bias_regularization` to `0.1`. We would have to add an entry to the config file such as:

```
anat_warp.bias_regularization: 0.1
```

To know the names of the nodes defined in `pypes`, ,
you will have to check either the **"Related settings"** sections in the 
description of the workflows here or the corresponding workflow source code.