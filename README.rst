Python bindings to TEOS-10
--------------------------

This package is still in beta testing phase and some of its functions are
likely to evolve over time.

Includes:

- Sea-Ice-Water (SIA) Library 3.0 + simple code fixes (for easy compilation)
- An extra module layer `Liq_Ice_Air_5a.F90`, `Liq_Ice_Air_5b.F90`, `Liq_Ice_Air_5c.F90` and
  `Liq_Ice_Air_5d.F90` with the temperature, the entropy, the potential temperature,
  the equivalent potential temperature and the pseudo equivalent potential temperature
  valid over the whole region of validity of TEOS-10.
- Two other extra modules `Liq_Air_4d.F90` and `Ice_Air_4d.F90` that are very similar
  to `Liq_Air_4b.F90` and `Ice_Air_4b.F90` but where one can specify a maximum
  saturation mass fraction in the definition of the Gibbs function. This allows
  to compute the equivalent potential temperature in a clean manner.
- A framework to add python bindings to the Fortran version of TEOS-10 with basic examples.

This package provides a simple way to import the TEOS-10 SIA Library into python.
The current beta version offers a very basic interface to the library. Future versions
will include a better interface and a more extensive physical description of the different
routines available.

NOTE:
    This package requires a Fortran 90 compiler. Sometimes numpy.disutils will not find
    the installed Fortran compiler and will return an error. It this happens, download the
    source code and run::
    
        $ f2py -c --help-fcompiler      

    To find the list of available compilers and then compile::

        $ python setup.py build --fcompiler=gnu95

    if 'gnu95' was in the list.

VERSION HISTORY
^^^^^^^^^^^^^^^

1.0.4 : New equivalent potential temperature and its approximate version.

1.0.3 : Must compute f2py in ``setup.py`` but now using a more portable method.

1.0.1 : Package comes with pre-computed f2py files for easier portability.

1.0 : Important bug fix on PYPI

0.9.5 : Some minimal testing and code publication improvements

0.9.2 : Finished moving netCDF interface to package ``pyteos_nc``

0.9.1 : Moved netCDF interface to package ``pyteos_nc``

0.9 : Improve netCDF interface

0.8 : Fixed a warning signal with missing values

0.7 : The interpolation is now much faster

0.6 : Faster netcdf interface for computing interpolated thermodynamic quantities.


Frederic Laliberte, February 17 2015
