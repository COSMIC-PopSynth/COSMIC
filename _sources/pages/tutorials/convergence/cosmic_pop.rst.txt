.. _fixedpop:

######################################
Generate a population the `COSMIC` way
######################################
Beyond providing a simple interface and several updates to BSE (`Hurley+2002 <https://ui.adsabs.harvard.edu/abs/2002MNRAS.329..897H/abstract>`_), COSMIC is designed to adapt the number of systems simulated to each binary evolution model a user selects. This is done by iteratively simulating binaries initialized with ZAMS binary parameters and a star formation history until the distributions of binary parameters specified by the user converge. This process is carried out by ``cosmic-pop``. Once this population is simulated, an astrophysical population can be sampled from the output of ``cosmic-pop``.


********************
cosmic-pop
********************
The `fixed population` is the output of ``cosmic-pop`` and is designed to fully describe the distribution of binary parameters that results from a user specified star formation history (SFH) and binary evolution model (BSEDict, specified in an inifile). The fixed population is only simulated once and only contains information about the intrinsic properties of the binary (e.g. masses, semimajor axes, metallicities, etc.) Information about the location and orientation of each binary is `not` contained in the fixed population. The arguments necessary to run the ``cosmic-pop`` executable can be found using the help command:

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

The ``rand_seed`` subsections allows you to initialize the population with a random seed for reproduceability. Note that for each simulated binary, COSMIC also returns the initial conditions, including a random seed that will reproduce that exact binary. The random seed produced in the rand_seed subsection allows full populations to be reproduced. This is particularly useful when comparing two popuations, e.g. binary black holes and binary neutron stars, which should be simulated separately, but using the same rand_seed value.

The ``sampling`` subsection allows you to control how you sample and achieve the initial binaries that you evolve.

The ``bse`` subsection is where all the BSE flags are specified.

-------------------
Sample command line
-------------------

Below is an example command line call for ``cosmic-pop`` with multiprocessing (i.e. shared memory on a single computer):

.. code-block:: bash

   cosmic-pop --final-kstar1 13 14 --final-kstar2 13 14 --inifile Params.ini --Nstep 1000 --Niter 1000000000 -n 2

Below is an example command line call for ``cosmic-pop`` with MPI (i.e. across many computers):

.. code-block:: bash

    mpiexec -n 2 python PATHTOEXECTUABLEDIRECTORY/cosmic-pop --final-kstar1 13 14 --final-kstar2 13 14 --inifile Params.ini --Nstep 1000 --Niter 1000000000

A breakdown of each argument follows:

* ``--final-kstar1 13 14 --final-kstar2 13 14`` : this tells COSMIC to keep all systems where the primary star is a BH or NS and the secondary star is also a BH or NS.

* ``--inifile Params.ini`` : this tells COSMIC where to look for the BSE flags that set the binary evolution model; in this case, the inifile is assumed to be in the same directory that the command line call is made

* ``--Nstep 5000 --Niter 10000000 -n 2`` : this tells COSMIC to use 2 processors (or in the MPI case, MPI tasks) to evolve a maximum of 1e7 systems and check in every 5000 systems to track how the shape of the distributions of the parameters specified in convergence-params change

===================
Stopping conditions
===================

There are two stopping conditions for ``cosmic-pop``:

1. The shape of the parameter distributions for the parameters specified in convergence_params section of the inifile file converge to a shape, regardless of adding new binaries to the evolved population. This is quantified by the match criteria detailed in Breivik+2020

2. The number of binaries sampled exceeds ``--Niter``.

==============================
Output of cosmic-pop
==============================

PLEASE SEE :ref:`output_info` for more information about the output data frames including
what each column means and the units.

The output of ``cosmic-pop`` is the `fixed population`, an hdf5 file with a naming scheme that tells you the Galactic component and final kstars of the population; the data file created by the ``cosmic-pop`` call above is: dat_DeltaBurst_13_14_13_14.h5.

The fixed population contains several pandas DataFrames accessed by the following keys:

* ``conv`` : The converged population whose parameters are sepcified by the ``convergence`` subsection of the inifile

* ``bpp`` : The evolutionary history of binaries which satisfy the user-specified final kstars and filter in the ``convergence`` subsection

* ``bcm`` : The final state of binaries in the bcm array which satisfy the user-specified final kstars and filter in the ``convergence`` subsection

* ``kick_info`` : The magnitude and direction of natal kicks, three dimensional systemic velocity changes, total tilt of orbital plane, and azimuthal angle of orbital angular momentum axis with respect to spins

* ``initC`` : The initial conditions for each binary which satisfies the user-specified final kstars and filter in the ``convergence`` subsection

* ``bpp_singles`` : The evolutionary history of single stars which satisfy the user-specified final kstars and filter in the ``convergence`` subsection

* ``bcm_singles`` : The final state of single stars in the bcm array which satisfy the user-specified final kstars and filter in the ``convergence`` subsection

* ``kick_info_singles`` : The magnitude and direction of natal kicks, three dimensional systemic velocity changes, total tilt of orbital plane, and azimuthal angle of orbital angular momentum axis with respect to spins

* ``initC_singles`` : The initial conditions for each single star which satisfies the user-specified final kstars and filter in the ``convergence`` subsection

* ``idx`` : An integer that keeps track of the total number of simulated binaries to maintain proper indexing across several runs of ``cosmic-pop``

* ``match`` : Tracks the convergence where match = Log :sub:`10` (1-convergence)

* ``mass_binaries`` : Tracks the total mass of binaries needed to create the fixed population

* ``mass_singles`` : Tracks the total mass of single stars needed to create the fixed population; if the binary fraction is 100%, the mass in singles will be zero

* ``mass_stars`` : Tracks the total mass of all stars, including binaries and singles, needed to create the fixed population

* ``n_binaries`` : Tracks the total number of binaries needed to create the fixed population

* ``n_singles`` : Tracks the total number of single stars needed to create the fixed population

* ``n_stars`` : Tracks the total number of stars, where n_stars = n_singles + 2*n_binaries, needed to create the fixed population

The ``conv``, ``bpp``, ``bcm``, and ``initCond`` DataFrames share a common column called ``bin_num`` which is used to index the population across the DataFrames.


**************************************
scaling to an astrophysical population
**************************************
Once the fixed population is simulated, you can scale the simulation to an astrophysical population by resampling the ``conv`` DataFrame. 

First, we need to load the data which is saved in the same directory where` ``cosmic-pop`` is called:

.. ipython::

    In [1]: import pandas

    In [2]: import numpy

    In [3]: conv = pandas.read_hdf('data/dat_DeltaBurst_13_14_13_14.h5', key='conv')

    In [4]: total_mass = pandas.read_hdf('data/dat_DeltaBurst_13_14_13_14.h5', key='mass_stars')

    In [5]: N_stars = pandas.read_hdf('data/dat_DeltaBurst_13_14_13_14.h5', key='n_stars')

.. note::

    The masses and numbers of stars/binaries is saved at each iteration, so you'll need to take the maximum mass/number:

.. ipython::

    In [6]: print(N_stars)

    In [7]: total_mass = max(numpy.array(total_mass))[0]

    In [8]: N_stars = max(numpy.array(N_stars))[0]

    In [9]: print(total_mass, N_stars)

Since COSMIC tracks both the total number and total mass of stars formed, you can either scale by number or mass. 

To scale by number, multiply the number of systems in the conv array by the ratio of the total number of stars in the astrophysical population to the total number of stars used to generate the population:

.. ipython::

    In [10]: N_astro = 1e10

    In [11]: N_13_14_13_14_astro = int(len(conv)*N_astro/N_stars)

    In [12]: print(N_13_14_13_14_astro)


Instead, if you want to scale by mass, you can choose between supplying your own mass or using the built in COSMIC models. The process for a user-supplied mass is nearly identical to scaling by number:

.. ipython::

    In [13]: M_astro = 1e11

    In [14]: N_13_14_13_14_astro = int(len(conv)*M_astro/total_mass)

    In [15]: print(N_13_14_13_14_astro)

Now we can generate the astrophysical population:

.. ipython::

    In [16]: pop_astro = conv.sample(N_13_14_13_14_astro, replace=True)

If you specified a star formation history in the inifile for the ``cosmic-pop`` call (i.e. by setting ``SF_start`` and ``SF_duration`` to follow a star formation history of your choice), the data in the resampled population will also be consistent with that star formation history. If you assigned the same ``SF_start`` to all binaries in the ``cosmic-pop`` call, you can assign birth times in post processing. As an example, we can assign a uniform birth time, assuming that the age of the population is 10 Gyr:

.. ipython::

    In [23]: pop_astro['tbirth'] = np.random.uniform(0, 10000, N_13_14_13_14_astro)

Since we are interested in NSs/BHs that form before the present, we can filter out anything that has a formation time after 10 Gyr: 

.. ipython::

    In [24]: pop_astro = pop_astro.loc[pop_astro.tbirth + pop_astro.tphys < 10000]

This leaves us with a population of NSs/BHs at formation where the formation time is the sum of the birth time and tphys in the ``conv`` array. 
