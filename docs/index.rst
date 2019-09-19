##################################
Welcome to COSMIC's documentation!
##################################
COSMIC is a rapid binary population synthesis suite with a special purpose of generating realistic compact binary populations. 

COSMIC currently implements binary evolutionary processes using BSE (`Hurley+2002 <http://adsabs.harvard.edu/abs/2002MNRAS.329..897H>`_). Several modifications have been applied to BSE to account for recent updates to binary evolution especially important to compact binary formation (e.g. metallicity-dependent stellar winds or black hole natal kick strengths). For a detailed discussion of these modifications, see Breivik+2019 in prep.

************
Using COSMIC
************
COSMIC's primary purpose is to generate synthetic populations with an adaptive size based on how the binary parameter disributions converge.
This is done through an executable that is installed when COSMIC is installed:

* cosmic-pop (see :ref:`fixedpop`)

COSMIC can also be used to simulate a single binary at a time, a list of multiple binaries, a grid of binaries, or a fixed population size. For more details on how to use COSMIC in these ways see :ref:`examples` and :ref:`runpop`.

*****************
Table of Contents
*****************

.. toctree::
   :maxdepth: 1

   install/index
   inifile/index
   output_info/index
   examples/index
   runpop/index
   fixedpop/index

*********************
Package documentation
*********************

Please consult these pages for more details on using COSMIC:

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
