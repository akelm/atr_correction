#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 19 10:23:18 2017

@author: ania
"""

from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import numpy # to get includes

# extensions = [Extension("*", ["yamsx/*.pyx"])]

setup(
    cmdclass = {'build_ext': build_ext},
    ext_modules = [Extension("kramers", ["kramers.pyx"], )],
    include_dirs = [numpy.get_include(),],
)