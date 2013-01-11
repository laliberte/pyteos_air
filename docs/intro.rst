============
Introduction
============
A quick overview of TEOS for wet air with ice and liquid
--------------------------------------------------------
This python package provides bindings to the Fortran version
of the `SIA Library version 3.0 <http://www.teos-10.org/>`_. 
It includes some code modifications for easy compiling and adds several modules to tackle the definition
of the Equivalent Potential Temperature.

This code has been developed on CentOS 5.6 on a x86_64 architecture using gcc 4.1.2.
It has been successfully compiled on CentOS 6.2 on a x86_64 architecture using ifort 12.1.3.
Finally, it was compiled on Mac OS X Lion using the ``gcc-lion.tar.gz`` (4.8) distribution from 
`HPC on Mac OS X <http://hpc.sourceforge.net/>`_.

Sometimes the Fortran 90 compiler is not found on a system. In this case, try::

    f2py -c --help-fcompiler

and use any one compiler in the list. For example, if ``g95`` is available compile with::

    python setup.py build --fcompiler=g95

Pickle
------
THe functions contained within this package can all be saved using pickle, making the use of the multiprocessing toolbox possible.

Frédéric B. Laliberté, November 20, 2012
