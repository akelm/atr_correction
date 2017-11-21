#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 19 10:23:18 2017

@author: ania
"""

from distutils.core import setup
from Cython.Build import cythonize

# extensions = [Extension("*", ["yamsx/*.pyx"])]

setup(
    ext_modules = cythonize("kramers.pyx")
)