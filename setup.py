#!/usr/bin/env python
"""Description:
Setup script for CRISPResso2 -- Software pipeline for the analysis of genome editing outcomes from deep sequencing data
@status:  beta
@version: $Revision$
@author:  Kendell Clement
@contact: kclement@mgh.harvard.edu

CRISPResso2 - Kendell Clement and Luca Pinello 2018
Software pipeline for the analysis of genome editing outcomes from deep sequencing data
(c) 2018 The General Hospital Corporation. All Rights Reserved.
"""

from setuptools import setup, Extension

# Use build_ext from Cython if found
command_classes = {}
try:
    import Cython.Distutils
    command_classes['build_ext'] = Cython.Distutils.build_ext
    has_cython = True
except:
    has_cython = False

ext = '.pyx' if has_cython else '.c'

from numpy import get_include as numpy_get_include
numpy_include_dir = [numpy_get_include()]

import sys
import re

def main():
    if float(sys.version[:3])<2.7 or float(sys.version[:3])>=2.8:
        sys.stderr.write("CRITICAL: Python version must be 2.7!\n")
        sys.exit(1)

    version = re.search(
    	'^__version__\s*=\s*"(.*)"',
    	open('CRISPResso2/CRISPRessoShared.py').read(),
    	re.M
    	).group(1)

    ext_modules = [
            Extension("CRISPResso2.CRISPRessoCOREResources", ["CRISPResso2/CRISPRessoCOREResources" + ext], include_dirs=numpy_include_dir, extra_compile_args=['-w','-Ofast'] ),
            Extension("CRISPResso2.CRISPResso2Align", ["CRISPResso2/CRISPResso2Align" + ext], include_dirs=numpy_include_dir, extra_compile_args=['-w','-Ofast'] ),
                       ]
    if has_cython:
        from Cython.Build import cythonize
        ext_modules = cythonize(ext_modules, language_level="2")

    setup(name="CRISPResso2",
          version=version,
          author='Kendell Clement',
          author_email='kclement@mgh.harvard.edu',
          url='http://github.com/pinellolab/CRISPResso2',
          package_dir={'CRISPResso2' : 'CRISPResso2'},
          include_package_data = True,
          packages=['CRISPResso2'],
      	  entry_points = {
        	"console_scripts": ['CRISPResso = CRISPResso2.CRISPRessoCORE:main',
          'CRISPRessoBatch = CRISPResso2.CRISPRessoBatchCORE:main',
          'CRISPRessoPooled = CRISPResso2.CRISPRessoPooledCORE:main',
          'CRISPRessoWGS = CRISPResso2.CRISPRessoWGSCORE:main',
          'CRISPRessoCompare = CRISPResso2.CRISPRessoCompareCORE:main',
          'CRISPRessoPooledWGSCompare = CRISPResso2.CRISPRessoPooledWGSCompareCORE:main',
              ]
           },
          description="Software pipeline for the analysis of genome editing outcomes from deep sequencing data",
          classifiers=[
              'Development Status :: 4 - Beta',
              'Environment :: Console',
              'Intended Audience :: Developers',
              'Intended Audience :: Science/Research',
              'License :: Other/Proprietary License',
              'Operating System :: POSIX',
              'Topic :: Scientific/Engineering :: Bio-Informatics',
              'Programming Language :: Python :: 2 :: Only',
              'Programming Language :: Cython',
              ],
          install_requires=[
              'numpy>=1.9',
              'pandas>=0.15',
              'matplotlib>=1.3.1',
              'biopython>=1.6.5',
              'argparse>=1.3',
			  'seaborn>=0.7.1',
              ],
          cmdclass = command_classes,
          ext_modules = ext_modules
          )

if __name__ == '__main__':
    main()
