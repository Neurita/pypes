#!/usr/bin/env python
"""
pypes
-----

Reusable neuroimaging pipelines using nipype
"""
import io

from setuptools import setup, find_packages

requirements = [
    'scipy>=1.0.0',
    'hansel>=2.0.1',
    'scikit_learn>=0.19.1',
    'matplotlib==2.1.0',
    'nilearn>=0.4.0',
    'kaptan==0.5.9',
    'nibabel>=2.2.1',
    'nipype>=1.0.2',
    'boyle>=0.2.0',
    'pandas>=0.22.0',
    'click>=6.7',
    'nipy>=0.4.2',
    'pydicom>=1.0.1',
]


# long description
def read(*filenames, **kwargs):
    encoding = kwargs.get('encoding', 'utf-8')
    sep = kwargs.get('sep', '\n')
    buf = []
    for filename in filenames:
        with io.open(filename, encoding=encoding) as f:
            buf.append(f.read())
    return sep.join(buf)


setup_dict = dict(
    name='neuro_pypes',
    version='1.1.2',
    description='Reusable and configurable neuroimaging pipelines with Nipype.',
    license='Apache License, Version 2.0',
    author='Alexandre Savio',
    author_email='alexsavio@gmail.com',
    maintainer='',
    maintainer_email='',
    packages=find_packages(exclude=['tests']),
    install_requires=requirements,
    scripts=[],
    entry_points='''
      [console_scripts]
      nitap=neuro_pypes.cli:cli
      ''',
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
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
    ],
    extras_require={
        'tests': [],
        'docs': ['mkdocs', 'recommonmark']
    }
)

if __name__ == '__main__':
    setup(**setup_dict)
