#!/usr/bin/env python

"""
Install the packages you have listed in the requirements file you input as
first argument.
"""

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from fabric.api import task, local, cd

import sys
import os
import os.path as op
import fileinput
import subprocess
import shutil

from   glob       import glob
from   setuptools import Command, setup, find_packages
from   pip.req    import parse_requirements

# Get version without importing, which avoids dependency issues
MODULE_NAME = find_packages(exclude=['tests'])[0]
VERSION_PYFILE = op.join(MODULE_NAME, 'version.py')
exec(compile(open(VERSION_PYFILE).read(), VERSION_PYFILE, 'exec'))

#get current dir
CWD = op.realpath(op.curdir)


def get_requirements(req_files):
    """Parse all requirements files given and return a list of the
       dependencies.

    Parameters
    ----------
    req_files: list of strings
        List of requirements file paths

    Returns
    -------
    dependencies: list of strings
        List of dependencies
    """
    install_deps = []

    if isinstance(req_files, str):
        req_files = [req_files]

    try:
        for fpath in args:
            install_deps.extend([str(d.req) for d in parse_requirements(fpath)])
    except:
        print('Error reading {} file looking for dependencies.'.format(fpath))

    return [dep for dep in install_deps if dep != 'None']


def recursive_glob(base_directory, regex=None):
    """
    Uses glob to find all files that match the regex in base_directory.

    Parameters
    ----------
    base_directory: str
        Root folder from where to start the search.

    regex: str
        Regular expression to be glob matched

    Returns
    -------
    files: set of str
        Set of filepaths that match the glob search.
    """
    if regex is None:
        regex = ''

    files = glob(os.path.join(base_directory, regex))
    for path, dirlist, filelist in os.walk(base_directory):
        for dir_name in dirlist:
            files.extend(glob(os.path.join(dir_name, regex)))

    return files


def recursive_remove(work_dir=CWD, regex='*'):
    [os.remove(fn) for fn in recursive_glob(work_dir, regex)]


def recursive_rmtrees(work_dir=CWD, regex='*'):
    [shutil.rmtree(fn, ignore_errors=True) for fn in recursive_glob(work_dir, regex)]


@task
def install_deps():
    deps = get_requirements(['requirements.txt'])

    try:
        for dep_name in deps:
            cmd = "pip install '{0}'".format(dep_name)
            print('#', cmd)
            subprocess.check_call(cmd, shell=True)
    except:
        print('Error installing {}'.format(dep_name))


@task
def version():
    print(__version__)


@task
def install():
    clean()
    install_deps()
    local('python setup.py install')


@task
def develop():
    clean()
    install_deps()
    local('python setup.py develop')


@task
def clean(work_dir=CWD):
    clean_build(work_dir)
    clean_pyc(work_dir)


@task
def clean_build(work_dir=CWD):
    shutil.rmtree('build', ignore_errors=True)
    shutil.rmtree('dist', ignore_errors=True)
    recursive_rmtrees(work_dir, '__pycache__')
    recursive_rmtrees(work_dir, '*.egg-info')


@task
def clean_pyc(work_dir=CWD):
    recursive_remove(work_dir, '*.pyc')
    recursive_remove(work_dir, '*.pyo')
    recursive_remove(work_dir, '*~')


@task
def lint():
    local('flake8 ' + MODULE_NAME + ' test')


@task
def test():
    local('py.test')


@task
def test_all():
    local('tox')


@task
def coverage():
    local('coverage local --source ' + MODULE_NAME + ' setup.py test')
    local('coverage report -m')
    local('coverage html')
    local('open htmlcov/index.html')


@task
def docs(doc_type='html'):
    os.remove(op.join('docs', MODULE_NAME + '.rst'))
    os.remove(op.join('docs', 'modules.rst'))
    local('sphinx-apidoc -o docs/ ' + MODULE_NAME)
    os.chdir('docs')
    local('make clean')
    local('make ' + doc_type)
    os.chdir(CWD)
    local('open docs/_build/html/index.html')


@task
def release():
    clean()
    local('python setup.py sdist upload')
    local('python setup.py bdist_wheel upload')


@task
def sdist():
    clean()
    local('python setup.py sdist')
    local('python setup.py bdist_wheel upload')
    print(os.listdir('dist'))
