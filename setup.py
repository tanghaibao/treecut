#!/usr/bin/env python
# -*- coding: UTF-8 -*-

from setuptools import setup
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


def import_init(filename="__init__.py"):
    """ Get various info from the package without importing them
    """
    import ast

    with open(filename) as init_file:
        module = ast.parse(init_file.read())

    itr = lambda x: (ast.literal_eval(node.value) for node in ast.walk(module) \
        if isinstance(node, ast.Assign) and node.targets[0].id == x)

    try:
        return next(itr("__author__")), \
               next(itr("__email__")), \
               next(itr("__license__")), \
               next(itr("__version__"))
    except StopIteration:
        raise ValueError("One of author, email, license, or version"
                    " cannot be found in {}".format(filename))


author, email, license, version = import_init(filename="treecut/__init__.py")

setup(
    name=name,
    version=version,
    author=author,
    author_email=email,
    packages=[name],
    scripts=glob('scripts/*.py') + ["treecut.py"],
    license=license,
    classifiers=classifiers,
    url='http://github.com/tanghaibao/treecut',
    description="Find nodes in hierarchical clustering that are statistically significant",
    long_description=open("README.md").read(),
    install_requires=['numpy', 'fisher', 'ete2'],
)
