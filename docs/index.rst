.. cosmic documentation master file, created by
   sphinx-quickstart on Thu Apr 21 14:05:08 2016.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

##################################
Welcome to cosmic's documentation!
##################################

*****************
Installing cosmic
*****************

The easiest method to install cosmic is using `pip <https://pip.pypa.io/en/stable/>`_ directly from the `GitHub repository <https://github.com/COSMIC-PopSynth/COSMIC.git>`_:

.. code-block:: bash

   $ pip install git+https://github.com/COSMIC-PopSynth/COSMIC.git


*****************
Table of Contents
*****************

.. toctree::
   :maxdepth: 4

   examples/index
   runpop/index

*****************
How to run cosmic
*****************

Here is an example commandline:

.. code-block:: bash
    runFixedPop --inifile Params.ini --Niter 1000000 --Nstep 100000 --galaxy-component ThinDisk --nproc=2 --final-kstar1=11 --final-kstar2=11

For more details see :ref:`command-line`.

*********************
Package documentation
*********************

Please consult these pages for more details on using cosmic:

.. toctree::
   :maxdepth: 1

   command-line/index

*******
Fortran
*******

evolv2
======

******************
Indices and tables
******************

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
