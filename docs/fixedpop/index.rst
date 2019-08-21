.. _fixedpop:

######################################
Generate a population the `cosmic` way
######################################
There are two executables that are installed when you pip install cosmic:

* cosmic-pop

* gxRealization

These two executables generate the two components of a cosmic Milky Way population.

********************
cosmic-pop
********************
The `fixed population` is simulated first and is designed to fully describe the population of binaries that results from a user specified star formation history (SFH) and binary evolution model (BSEDict, specified in an inifile). The fixed population is only simulated once and only contains information about the intrinsic properties of the binary (e.g. masses, semimajor axes, metallicities, etc.) Information about the location in the Galaxy of each binary is `not` contained in the fixed population. The arguments necessary to run the cosmic-pop executable can be found using the help command:

.. code-block:: bash

   cosmic-pop -h

.. program-output:: cosmic-pop -h

======
Inputs
======

----------
Params.ini
----------

PLEASE SEE :ref:`inifile` for detailed information about the inifile and how it is constructed

The inifile contains five subsections: ``filters``, ``convergence``, ``rand_seed``, and ``bse``.

The ``filters`` subsection allows you to specify how you would like to filter the binary population. See the inifile for a description of each flag.

The ``convergence`` subsection allows you to specify a particular region of parameter space where you would like the convergence of each binary parameter distribution to be tested. The only implemented convergence filter is for LISA binaries, where the convergence is computed only for binaries with orbital period less than 5000 s. This allows for low probability, but high signal to noise binaries with very short orbital periods to be fully accounted for.

The ``rand_seed`` subsections allows you to initialize the population with a random seed for reproduceability. Note that for each simulated binary, cosmic also returns the initial conditions, including a random seed that will reproduce that exact binary. The random seed produced in the rand_seed subsection allows full populations to be reproduced. This is particularly useful when comparing two popuations, e.g. binary black holes and binary neutron stars, which should be simulated separately, but using the same rand_seed value.

The ``sampling`` subsection allows you to control how you sample and achieve the initial binaries that you evolve.

The ``bse`` subsection is where all the bse flags are specified.

-------------------
Sample command line
-------------------

Below is an example command line call for cosmic-pop:

.. code-block:: bash

   cosmic-pop --final-kstar1 13 14 --final-kstar2 13 14 --inifile Params.ini --Nstep 1000 --Niter 1000000000 -n 2

A breakdown of each argument follows:

* ``--final-kstar1 13 14 --final-kstar2 13 14`` : this tells cosmic to keep all systems where the primary star is a BH or NS and the secondary star is also a BH or NS.

* ``--inifile Params.ini`` : this tells cosmic where to look for the BSE flags that set the binary evolution model; in this case, the inifile is assumed to be in the same directory that the command line call is made

* ``--Nstep 5000 --Niter 10000000 -n 1`` : this tells cosmic to use 1 processor to evolve a maximum of 1e7 systems and check in every 5000 systems to track how the shape of the distributions of the parameters specified in convergence-params change

===================
Stopping conditions
===================

There are two stopping conditions for ``cosmic-pop``:

1. The shape of the parameter distributions for the parameters specified in convergence_params section of the inifile file converge to a shape, regardless of adding new binaries to the evolved population. This is quantified by the match criteria detailed in Breivik+2019 in prep

2. The number of binaries sampled exceeds ``--Niter``.

==============================
Output of cosmic-pop
==============================

PLEASE SEE :ref:`output_info` for more information about the output data frames including
what each column means and the units.

The output of ``cosmic-pop`` is the `fixed population`, an hdf5 file with a naming scheme that tells you the Galactic component and final kstars of the population; the data file created by the cosmic-pop call above is: dat_DeltaBurst_13_14_13_14.h5.

The fixed population contains three pandas DataFrames accessed by the following keys:

* ``bcm`` : The final state of the converged population at the present epoch

* ``bpp`` : The evolutionary history of the systems in the bcm data set at key moments in the evolution

* ``initCond`` : The initial conditions for each binary in the bcm data set

Each of these DataFrames shares a common column called ``bin_num`` which is used to index the population across the DataFrames.
