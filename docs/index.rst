##################################
Welcome to cosmic's documentation!
##################################
cosmic is a rapid binary population synthesis suite with a special purpose of generating realistic Milky Way compact binary populations. 

cosmic currently implements binary evolutionary processes using BSE (`Hurley+2002 <http://adsabs.harvard.edu/abs/2002MNRAS.329..897H>`_). Several modifications have been applied to BSE to account for recent updates to binary evolution especially important to compact binary formation (e.g. metallicity-dependent stellar winds or black hole natal kick strengths). For a detailed discussion of these modifications, see Breivik+2018 in prep.

*****************
Installing cosmic
*****************

The easiest method to install cosmic is using `pip <https://pip.pypa.io/en/stable/>`_ directly from the `GitHub repository <https://github.com/COSMIC-PopSynth/COSMIC.git>`_:

-----
MacOS
-----

.. code-block:: bash

    conda create --name cosmic python=3.6
    source activate cosmic
    conda install gcc numpy
    conda install -c uvcdat gfortran
    pip install git+https://github.com/COSMIC-PopSynth/COSMIC.git

----
Unix
----

.. code-block:: bash

    conda create --name cosmic python=3.6
    source activate cosmic
    conda install numpy
    pip install git+https://github.com/COSMIC-PopSynth/COSMIC.git


*******************************
Using cosmic to evolve binaries
*******************************

For more details on how to use cosmic, see :ref:`examples`, :ref:`runpop`, and :ref:`fixedpop`.

How to run a fixed population
-----------------------------

One of the main product of this package is the command-line executable `runFixedPop`,
which takes an excess noise time makes an omega scan of the event and classifies the image.

To run an analysis:

.. code-block:: bash

   $ runFixedPop --final_kstar1 11 --final_kstar2 11 --galaxy_component ThinDisk --initial_samp multidim --Nstep 500 --Niter 10000 -n 1 --inifile my-config-file.ini

where ``./my-config-file.ini`` is the path of your BSE configuration file.
In the `examples` folder is an example of an ini file.

For a full list of command-line argument and options, run

.. code-block:: bash

   $ runFixedPop --help

For more details see :ref:`runFixedPop`.


*****************
Table of Contents
*****************

.. toctree::
   :maxdepth: 2

   install/index
   examples/index
   runpop/index
   fixedpop/index

*********************
Package documentation
*********************

Please consult these pages for more details on using cosmic:

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
