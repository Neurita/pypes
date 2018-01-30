#!/usr/bin/env python
"""
pypes
-----

Reusable neuroimaging pipelines using nipype
"""

from __future__ import print_function

import os
import io

from setuptools import setup, find_packages


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
MODULE_NAME = find_packages(exclude=['*.tests', '*.test.*'])[0]
VERSION_PYFILE = os.path.join(MODULE_NAME, 'version.py')
# set __version__ variable
exec(compile(read(VERSION_PYFILE), VERSION_PYFILE, 'exec'))

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
    install_requires=[
        # 'numpy>=1.13.3',
        'scipy>=1.0.0',
        'hansel>=0.9.5',
        'scikit_learn>=0.19.1',
        'matplotlib>=2.1.0',
        'nilearn>=0.3.1',
        'kaptan>=0.5.9',
        'nibabel>=2.2.1',
        'nipype>=0.14.0',
        'boyle>=0.2.0',
        'pandas>=0.21.0',
        'nipy',
        'pydicom',
    ],
    dependency_links=[
        'git+https://github.com/darcymason/pydicom@ebf6a79602348d003a1d1324c66626f9f2b05432#egg=pydicom',
        'git+https://github.com/nipy/nipy.git@d49e8292adad6619e3dac710752131b567efe90e#egg=nipy'
    ],
    scripts=[],
    long_description=read('README.md', 'CHANGES.md'),
    platforms='Linux/MacOSX',
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
    extras_require={
        'tests': [],
        'docs': ['mkdocs', 'recommonmark']
    })

if __name__ == '__main__':
    setup(**setup_dict)
