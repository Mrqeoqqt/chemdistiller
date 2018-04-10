# -*- coding: utf-8 -*-
"""
Created on Mon Aug 22 23:35:01 2016

@author: ilaponog
"""

from distutils.core import setup
from Cython.Build import cythonize

setup(
    ext_modules = cythonize("ProcessFragmentFile_V1.pyx")
)