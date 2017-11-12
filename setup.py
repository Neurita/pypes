#!/usr/bin/env python

"""
pypes
-----

Reusable neuroimaging pipelines using nipype
"""

from __future__ import print_function

import os.path as op
import io
import sys

from   setuptools              import setup, find_packages
from   setuptools.command.test import test as TestCommand


# long description
def read(*filenames, **kwargs):
    encoding = kwargs.get('encoding', 'utf-8')
    sep = kwargs.get('sep', '\n')
    buf = []
    for filename in filenames:
        with io.open(filename, encoding=encoding) as f:
            buf.append(f.read())
    return sep.join(buf)


# Get version without importing, which avoids dependency issues
MODULE_NAME = find_packages(exclude=['tests'])[0]
VERSION_PYFILE = op.join(MODULE_NAME, 'version.py')
# set __version__ variable
exec(compile(read(VERSION_PYFILE), VERSION_PYFILE, 'exec'))

dependencies = read('requirements.txt').split()
requires = [dep for dep in dependencies if not dep.startswith(('git', 'http'))]
links = [dep for dep in dependencies if dep.startswith(('git', 'http'))]
requires += [req.split('#egg=')[1] for req in links]

LICENSE = 'Apache License, Version 2.0'


setup_dict = dict(
    name=MODULE_NAME,
    version=__version__,
    description='Reusable and configurable neuroimaging pipelines with Nipype.',

    license=LICENSE,
    author='Alexandre Savio',
    author_email='alexsavio@gmail.com',
    maintainer='',
    maintainer_email='',

    packages=find_packages(),

    install_requires=requires,

    dependency_links=links,

    scripts=[],

    long_description=read('README.md', 'CHANGES.md'),

    platforms='Linux/MacOSX',

    # https://pypi.python.org/pypi?%3Aaction=list_classifiers
    classifiers=[
        'Programming Language :: Python',
        'Development Status :: 4 - Beta',
        'Natural Language :: English',
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: Apache Software License',
        'Operating System :: OS Independent',
        'Topic :: Scientific/Engineering :: Information Analysis',
        'Topic :: Scientific/Engineering :: Medical Science Apps.',
        'Operating System :: POSIX',
        'Operating System :: Unix',
        'Operating System :: MacOS',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
    ],

    extras_require={'tests': [],
                    'docs': ['mkdocs',
                             'recommonmark']
                   }
)

# Python3 support keywords
if sys.version_info >= (3,):
    setup_dict['use_2to3'] = False
    setup_dict['convert_2to3_doctests'] = ['']
    setup_dict['use_2to3_fixers'] = ['']


class PyTest(TestCommand):
    user_options = [('pytest-args=', 'a', "Arguments to pass to py.test")]

    def initialize_options(self):
        TestCommand.initialize_options(self)
        #self.pytest_args = []

    def finalize_options(self):
        TestCommand.finalize_options(self)
        self.test_suite = True

    def run_tests(self):
        # import here, cause outside the eggs aren't loaded
        import pytest
        errno = pytest.main(self.pytest_args)
        sys.exit(errno)

if __name__ == '__main__':
    setup(**setup_dict)
