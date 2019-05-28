#!/usr/bin/env python
"""
pypes
-----

Reusable neuroimaging pipelines using nipype
"""
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


requires, extras_requires = [], []

setup_dict = dict(
    name='neuro_pypes',
    version='1.2.0',
    description='Reusable and configurable neuroimaging pipelines with Nipype.',
    license='Apache License, Version 2.0',
    author='Alexandre Savio',
    author_email='alexsavio@gmail.com',
    maintainer='',
    maintainer_email='',
    packages=find_packages(exclude=['tests']),
    install_requires=requires,
    extra_requires=extras_requires,
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
