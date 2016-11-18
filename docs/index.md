# Welcome to Pypes

Pypes is a Python module for brain PET and multimodal MRI processing.
It contains specific pipelines to pre-process this type of images.

The objectives of this module are:

- easy-to-use, complete, reusable, and configurable pre-processing pipelines, and
- high-quality results with the minimum manual intervention in the images.

First you will need to convert those DICOM files to NifTI files.
I recommend you to use [dcm2niix](https://github.com/rordenlab/dcm2niix) for this.
The you will probably want to clean and rename those files to have a homogeneous database.

I've tested these on ~140 subjects with very good results.
Without modifying parameters, 95% of the images were correctly registered and processed.
Very few subjects needed manual editings such as neck removal or AC-PC reorientation.

The spatial normalization here are only done with [SPM12](http://www.fil.ion.ucl.ac.uk/spm/software/spm12/) for now.


## Installation

The first thing you need to install are the external dependencies.

### External dependencies

Each pipeline type has a sub-set of external dependencies.
To install these, please check how to install them in their own manuals.
You can also find very helpful instructions [here](http://miykael.github.io/nipype-beginner-s-guide/installation.html).
If you are familiar with Bash, you can check our [neuro-docker](https://github.com/Neurita/neuro_docker).
Following, for each modality of MRI we list the dependencies you will need to have installed.

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
