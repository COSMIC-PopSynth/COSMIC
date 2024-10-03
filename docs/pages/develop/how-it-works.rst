*****************
Codebase overview
*****************

COSMIC is based on BSE and so is written primarily in Fortran. This underlying code is then accessed through
Python bindings, such that the end user has no need to interact directly with the Fortran code but still
benefits from the speed of the compiled code.

**Fortran code** -
The source code for the Fortran is found in the ``src/cosmic/src`` directory. This directory follows a similar
structure to the original BSE code, with separate files that handle different aspects of binary evolution. Of
particular note is ``evolv2.f`` which contains the main binary evolution code.

**Python binding** -
The Python bindings are found in the
``src/cosmic`` directory and are designed using `f2py <https://numpy.org/doc/stable/f2py/>`_. The functions
in the Fortran code are added to a shared object library which is then imported into Python with
``from cosmic import _evolvebin``.

**Sampling** -
There are classes and functions for sampling populations of sources in the ``src/cosmic/sample`` directory, which
is all written in Python.