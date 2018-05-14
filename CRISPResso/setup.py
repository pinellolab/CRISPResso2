'''
CRISPResso2 - Kendell Clement and Luca Pinello 2018
Software pipeline for the analysis of genome editing outcomes from deep sequencing data
(c) 2017 The General Hospital Corporation. All Rights Reserved.
'''
from distutils.core import setup
from Cython.Build import cythonize
import numpy

setup(
    ext_modules=cythonize("*.pyx"),
    include_dirs=[numpy.get_include()]
)
