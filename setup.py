"""Cython extension build configuration for CRISPResso2.

All package metadata, dependencies, and entry points are defined in
pyproject.toml. This file only exists to configure the Cython/C
extension modules, which setuptools does not support declaratively.
"""

from setuptools import setup, Extension, find_packages
from io import open

try:
    from Cython.Build import cythonize
    from Cython.Distutils import build_ext

    ext = '.pyx'
    has_cython = True
except ImportError:
    ext = '.c'
    has_cython = False

from numpy import get_include as numpy_get_include

ext_modules = [
    Extension(
        "CRISPResso2.CRISPRessoCOREResources",
        ["CRISPResso2/CRISPRessoCOREResources" + ext],
        include_dirs=[numpy_get_include()],
        extra_compile_args=['-w', '-Ofast'],
    ),
    Extension(
        "CRISPResso2.CRISPResso2Align",
        ["CRISPResso2/CRISPResso2Align" + ext],
        include_dirs=[numpy_get_include()],
        extra_compile_args=['-w', '-Ofast'],
    ),
]

if has_cython:
    ext_modules = cythonize(ext_modules, language_level="3")

setup(
    ext_modules=ext_modules,
    cmdclass={'build_ext': build_ext} if has_cython else {},
)
