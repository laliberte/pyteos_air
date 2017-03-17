# This Python file uses the following encoding: utf-8
import os
import sys
from setuptools import setup, find_packages
from numpy.distutils.core import setup, Extension

# Utility function to read the README file.
# Used for the long_description.  It's nice, because now 1) we have a top level
# README file and 2) it's easier to type in the README file than to put a raw
# string in below ...
def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

teos_files = [
'Constants_0.F90',
'Convert_0.F90',
'Maths_0.F90',
'Flu_1.F90',
'Ice_1.F90',
'Sal_1.F90',
'Air_1.F90',
'Flu_2.F90',
'Ice_2.F90',
'Sal_2.F90',
'Air_2.F90',
'Flu_3a.F90',
'Flu_3b.F90',
'Sea_3a.F90',
'Sea_3b.F90',
'Sea_3c.F90',
'Sea_3d.F90',
'Air_3a.F90',
'Air_3b.F90',
'Air_3c.F90',
'Liq_Vap_4.F90',
'Ice_Vap_4.F90',
'Ice_Liq_4.F90',
'Liq_Ice_Air_4.F90',
'Sea_Vap_4.F90',
'Sea_Liq_4.F90',
'Sea_Ice_4.F90',
'Sea_Ice_Vap_4.F90',
'Liq_Air_4a.F90',
'Sea_Air_4.F90',
'Liq_Air_4b.F90',
'Liq_Air_4c.F90',
'Liq_Air_4d.F90',
'Ice_Air_4a.F90',
'Ice_Air_4b.F90',
'Ice_Air_4c.F90',
'Ice_Air_4d.F90',
#'Liq_Ice_Air_5.F90',
'Liq_Ice_Air_5a.F90',
'Liq_Ice_Air_5b.F90',
'Liq_Ice_Air_5c.F90',
'Liq_Ice_Air_5d.F90',
#'Liq_Ice_Air_5e.F90',
'Flu_IF97_5.F90',
'Ice_Flu_5.F90',
'Liq_F03_5.F90',
'Sea_5a.F90',
'GSW_Library_5.F90',
'Convert_5.F90',
'Air_5.F90',
'OS2008_5.F90']
teos_files = ['src/SIA_library/' + file for file in teos_files]

sign_file = 'src/teos_air.pyf'
sources_interface = ['pyteos_air/liq_ice_air/liq_ice_air_f2py.F90']
teos_files.extend(sources_interface)

#BUILD TEOS-10. IT REQUIRES A SIGNATURE FILE THAT CAN BE GENERATED BY:

# generate signature file if it does not exist:
if not os.path.exists(sign_file):
    os.system('f2py -m teos_air --overwrite-signature -h ' + sign_file + ' ' +
              ' '.join(sources_interface))
teos_files.append(sign_file)

package_name = 'pyteos_air'
packages = [package_name + '.' + name for name in find_packages(exclude=['test'])]
packages.append(package_name)

setup(name=package_name,
      version="1.0.2",
      packages=find_packages(exclude=['test']),
      # metadata for upload to PyPI
      author=("F. Laliberte, R. Feistel, D. G. Wright, D. R. Jackett, K. Miyagawa, "
              "J. H. Reissmann, W. Wagner, U. Overhoff, C. Guder, A. Feistel, G. M. Marion"),
      author_email="frederic.laliberte@utoronto.ca",
      description="This package provides python bindings to the atmospheric component of TEOS-10. \
                   Also includes expanded metorological quantities like the equivalent potential temperature.",
      license="BSD",
      keywords="thermodynamics atmosphere climate",
      url="",   # project home page, if any
      classifiers=[
            "Development Status :: 4 - Beta",
            "Intended Audience :: Science/Research",
            "Natural Language :: English",
            "License :: OSI Approved :: BSD License",
            "Programming Language :: Python :: 2.7",
            "Programming Language :: Python :: 3.4",
            "Programming Language :: Python :: 3.5",
            "Programming Language :: Python :: 3.6",
            "Programming Language :: Fortran",
            "Topic :: Scientific/Engineering :: Atmospheric Science",
            "Topic :: Scientific/Engineering :: Mathematics"],
      long_description=read('README.rst'),
      install_requires=['numpy>=1.6', 'scipy'],
      extras_require={'testing': ['flake8',
                                  'coverage',
                                  'pytest-cov',
                                  'pytest']},
      zip_safe=False,
      ext_package='pyteos_air.teos_air',
      ext_modules=[Extension('teos_air', teos_files)])
