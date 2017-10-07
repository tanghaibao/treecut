#!/usr/bin/env python
# -*- coding: UTF-8 -*-

from setuptools import setup
from setup_helper import SetupHelper
from glob import glob


name = 'treecut'
classifiers = [
    'Development Status :: 4 - Beta',
    'Intended Audience :: Science/Research',
    'License :: OSI Approved :: BSD License',
    'Programming Language :: Python',
    'Programming Language :: Python :: 2',
    'Topic :: Scientific/Engineering :: Bio-Informatics',
    ]
with open('requirements.txt') as f:
    required = f.read().splitlines()

# Use the helper
h = SetupHelper(initfile="treecut/__init__.py", readmefile="README.md")
h.check_version(majorv=2, minorv=7)

setup(
    name=name,
    version=h.version,
    author=h.author,
    author_email=h.email,
    license=h.license,
    long_description=h.long_description,
    packages=[name],
    scripts=glob('scripts/*.py') + ["treecut.py"],
    classifiers=classifiers,
    url='http://github.com/tanghaibao/treecut',
    description="Find nodes in hierarchical clustering that are statistically significant",
    install_requires=required + ["ete2"],
)
