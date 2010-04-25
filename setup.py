#!/usr/bin/env python
# -*- coding: UTF-8 -*-

from distutils.core import setup #, Extension


setup(
      name="treecut",
      packages=['treecut'],
      requires=['fisher', 'ete2', 'statlib'],
      )
