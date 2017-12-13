#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (C) Katie Breivik (2017)
#
# This file is part of the aCOSMIC python package.
#
# aCOSMIC is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# aCOSMIC is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with aCOSMIC.  If not, see <http://www.gnu.org/licenses/>.

"""Setup the aCOSMIC package
"""

from __future__ import print_function

import sys
if sys.version < '2.6':
    raise ImportError("Python versions older than 2.6 are not supported.")

import glob
import os.path

from setuptools import find_packages

from numpy.distutils.core import setup, Extension

# set basic metadata
PACKAGENAME = 'aCOSMIC'
DISTNAME = 'aCOSMIC'
AUTHOR = 'Katie Breivik'
AUTHOR_EMAIL = 'katie.breivik@gmail.com'
LICENSE = 'GPLv3'

cmdclass = {}

# -- versioning ---------------------------------------------------------------

import versioneer
__version__ = versioneer.get_version()
cmdclass.update(versioneer.get_cmdclass())

# -- documentation ------------------------------------------------------------

# import sphinx commands
try:
    from sphinx.setup_command import BuildDoc
except ImportError:
    pass
else:
    cmdclass['build_sphinx'] = BuildDoc

# -- dependencies -------------------------------------------------------------

setup_requires = [
    'setuptools',
    'pytest-runner',
]
install_requires = [
    'numpy',
    'scipy',
    'astropy',
    'matplotlib',
    'gwpy',
    'pandas',
]
tests_require = [
    'pytest'
]
if sys.version_info < (2, 7):
    tests_require.append('unittest2')
extras_require = {
    'doc': [
        'sphinx',
        'sphinx-fortran',
        'numpydoc',
        'sphinx_rtd_theme',
        'sphinxcontrib_programoutput',
        'sphinxcontrib_epydoc',
    ],
}
# fortran compile

wrapper = Extension('aCOSMIC._evolvebin', sources=['aCOSMIC/src/comenv.f', 'aCOSMIC/src/corerd.f', 'aCOSMIC/src/deltat.f', 'aCOSMIC/src/dgcore.f', 'aCOSMIC/src/evolv2.f', 'aCOSMIC/src/gntage.f', 'aCOSMIC/src/hrdiag.f', 'aCOSMIC/src/instar.f', 'aCOSMIC/src/kick.f', 'aCOSMIC/src/mix.f', 'aCOSMIC/src/mlwind.f', 'aCOSMIC/src/mrenv.f', 'aCOSMIC/src/ran3.f', 'aCOSMIC/src/rl.f', 'aCOSMIC/src/star.f', 'aCOSMIC/src/zcnsts.f', 'aCOSMIC/src/zfuncs.f'], extra_compile_args = ["-O -g"])


# -- run setup ----------------------------------------------------------------

packagenames = find_packages()
scripts = glob.glob(os.path.join('bin', '*')) + glob.glob('data/*')

setup(name=DISTNAME,
      provides=[PACKAGENAME],
      version=__version__,
      description=None,
      long_description=None,
      ext_modules = [wrapper],
      author=AUTHOR,
      author_email=AUTHOR_EMAIL,
      license=LICENSE,
      packages=packagenames,
      include_package_data=True,
      cmdclass=cmdclass,
      scripts=scripts,
      setup_requires=setup_requires,
      install_requires=install_requires,
      tests_require=tests_require,
      extras_require=extras_require,
      use_2to3=True,
      classifiers=[
          'Programming Language :: Python',
          'Development Status :: 3 - Alpha',
          'Intended Audience :: Science/Research',
          'Intended Audience :: End Users/Desktop',
          'Intended Audience :: Science/Research',
          'Natural Language :: English',
          'Topic :: Scientific/Engineering',
          'Topic :: Scientific/Engineering :: Astronomy',
          'Topic :: Scientific/Engineering :: Physics',
          'Operating System :: POSIX',
          'Operating System :: Unix',
          'Operating System :: MacOS',
          'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
      ],
)
