#!/usr/bin/env python
# -*- coding: UTF-8 -*-

from setuptools import setup
from glob import glob

classifiers = [
    'Development Status :: 4 - Beta',
    'Intended Audience :: Science/Research',
    'License :: OSI Approved :: BSD License',
    'Programming Language :: Python',
    'Programming Language :: Python :: 2',
    'Topic :: Scientific/Engineering :: Bio-Informatics',
    ]

exec(open("treecut/version.py").read())
setup(
    name="treecut",
    version=__version__,
    author='Haibao Tang',
    author_email='tanghaibao@gmail.com',
    packages=['treecut'],
    scripts=glob('scripts/*.py'),
    license='BSD',
    classifiers=classifiers,
    url='http://github.com/tanghaibao/treecut',
    description="Find nodes in hierarchical clustering that are statistically significant",
    long_description=open("README.rst").read(),
    install_requires=['fisher', 'ete2', 'statlib'],
    )
