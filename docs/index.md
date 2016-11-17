# Welcome to Pypes


## Installation

The first thing you need to install are the external dependencies.

### External dependencies

Each pipeline type has a sub-set of external dependencies. To install these, please check how to install them in their own manuals.
You can also find very helpful instructions [here](http://miykael.github.io/nipype-beginner-s-guide/installation.html).

#### Anatomical MRI

-   [ANTs](http://stnava.github.io/ANTs/) (N4 inhomogeneity correction algorithm)
-   [SPM12](http://www.fil.ion.ucl.ac.uk/spm/software/spm12/) (anatomical image and atlas warping, and tissue segmentation).

#### FDG-PET

-   [SPM12](http://www.fil.ion.ucl.ac.uk/spm/software/spm12/) (anatomical co-registration and atlas normalization) and
-   [PETPVC](https://github.com/UCL/PETPVC) (PET partial volume correction).

#### DTI and tractography

-   [SPM12](http://www.fil.ion.ucl.ac.uk/spm/software/spm12/) (anatomical co-registration and atlas normalization),
-   [FSL](http://fsl.fmrib.ox.ac.uk/fsl/fslwiki/) ([Eddy](http://fsl.fmrib.ox.ac.uk/fsl/fslwiki/eddy) for motion and eddy-currents correction and fslmaths) and
-   [Camino](http://camino.cs.ucl.ac.uk/) (diffusion tensor model fitting and deterministic tractography).

Python dependency (installed later):

-   [Dipy](http://nipy.org/dipy/) (NL-means)

#### Resting-state fMRI preprocessing

-   [SPM12](http://www.fil.ion.ucl.ac.uk/spm/software/spm12/) (anatomical co-registration and atlas normalization),
-   [FSL](http://fsl.fmrib.ox.ac.uk/fsl/fslwiki/) (GLM for nuisance correction and fslmaths).

Python dependency (installed later):

-   [NiPy](http://nipy.org/nipy/documentation.html) (Realign for motion correction, this is already in the requirements.txt file) and

### Python dependencies

This package uses setuptools. First you may need to install the Python dependencies:

```shell
pip install -r requirements.txt
```

Then you can install it running:

```shell
python setup.py install
```

If you already have the dependencies listed in requirements.txt installed, to install in your home directory, use:

```shell
python setup.py install --user
```

To install for all users on Unix/Linux:

```shell
python setup.py build
sudo python setup.py install
```

You can also install it in development mode with:

```shell
python setup.py develop
```

## [Find me on Github](https://github.com/Neurita/pypes)

You can check the latest sources with the command:

```shell
git clone "https://www.github.com/neurita/pypes.git"
```

or if you have write privileges:

```shell
git clone git@github.com:neurita/pypes.git
```

If you are going to create patches for this project, create a branch for it from the master branch.

The stable releases are tagged in the repository.


## License

This software is licensed under the [Apache License, Version 2.0](http://www.apache.org/licenses/LICENSE-2.0).
