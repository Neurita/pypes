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

Please see the requirements.txt and pip_requirements.txt file.

Install
=======

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
