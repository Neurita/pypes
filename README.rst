.. -*- mode: rst -*-

pypes
=====

Reusable neuroimaging pipelines using nipype

.. image:: https://secure.travis-ci.org/neurita/pypes.png?branch=master
    :target: https://travis-ci.org/neurita/pypes

.. image:: https://coveralls.io/repos/neurita/pypes/badge.png
    :target: https://coveralls.io/r/neurita/pypes


Dependencies
============

Please see the requirements.txt file.

Apart from the requirements, the external dependencies for each pipeline are:


Anatomical MRI
--------------

- `SPM12 <http://www.fil.ion.ucl.ac.uk/spm/software/spm12/>`_ (anatomical image and atlas warping, and tissue segmentation).


FDG-PET
-------

- `SPM12 <http://www.fil.ion.ucl.ac.uk/spm/software/spm12/>`_ (anatomical co-registration and atlas normalization) and
- `PETPVC <https://github.com/UCL/PETPVC>`_ (partial volume correction).


DTI and tractography
--------------------

- `SPM12 <http://www.fil.ion.ucl.ac.uk/spm/software/spm12/>`_ (anatomical co-registration and atlas normalization),
- `FSL <http://fsl.fmrib.ox.ac.uk/fsl/fslwiki/>`_ (`Eddy <http://fsl.fmrib.ox.ac.uk/fsl/fslwiki/eddy>`_ for motion and eddy-currents correction and fslmaths) and
- `Camino <http://camino.cs.ucl.ac.uk/>`_ (diffusion tensor model fitting and deterministic tractography).

Resting-state fMRI preprocessing
--------------------------------

- `SPM12 <http://www.fil.ion.ucl.ac.uk/spm/software/spm12/>`_ (anatomical co-registration and atlas normalization),
- `NiPy <http://nipy.org/nipy/documentation.html>`_ (Realign for motion correction, this is already in the requirements.txt file) and
- `FSL <http://fsl.fmrib.ox.ac.uk/fsl/fslwiki/>`_ (GLM for nuisance correction and fslmaths).


Install
=======

To install the external dependencies, please check how to install them in their own manuals.


This package uses setuptools. You can install it running:

    python setup.py install

If you hve problems during this installation. First you may need to install the dependencies:

    pip install -r requirements.txt

If you already have the dependencies listed in requirements.txt installed,
to install in your home directory, use::

    python setup.py install --user

To install for all users on Unix/Linux::

    python setup.py build
    sudo python setup.py install

You can also install it in development mode with::

    python setup.py develop


Development
===========

Code
----

Github
~~~~~~

You can check the latest sources with the command::

    git clone https://www.github.com/neurita/pypes.git

or if you have write privileges::

    git clone git@github.com:neurita/pypes.git

If you are going to create patches for this project, create a branch for it
from the master branch.

The stable releases are tagged in the repository.


Testing
-------

TODO
