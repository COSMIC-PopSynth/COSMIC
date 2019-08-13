.. _fixedpop:

######################################
Generate a population the `cosmic` way
######################################
There are two executables that are installed when you pip install cosmic:

* runFixedPop

* gxRealization

These two executables generate the two components of a cosmic Milky Way population. 

***********
runFixedPop
***********
The `fixed population` is simulated first and is designed to fully describe the population of binaries that results from a user specified star formation history (SFH) and binary evolution model (BSEDict, specified in an inifile). The fixed population is only simulated once and only contains information about the intrinsic properties of the binary (e.g. masses, semimajor axes, metallicities, etc.) Information about the location in the Galaxy of each binary is `not` contained in the fixed population. The arguments necessary to run the runFixedPop executable can be found using the help command:

.. code-block:: bash

   runFixedPop -h 

.. program-output:: runFixedPop -h 

======
Inputs
======

------------------------
Star Formation Histories
------------------------

There are currently three SFHs implemented in cosmic, one each for the Milky Way thin disk, thick disk and bulge. Their assumptions are as follows:

* ThinDisk : constant star formation over 10 Gyr; the metallicity for the thin disk is traditionally set to solar

* ThickDisk : a 1 Gyr burst of constant star formation 11 Gyr in the past; the metallicity for the thick disk is traditionally set to be subsolar, we reccommend 15% solar

* Bulge : a 1 Gyr burst of constant 10 Gyr in the past; the metallicity for the bulge is traditionally set to solar 

New implementations of SFHs will be added to the sample_SFH modules in the independent and multidim samplers.  

----------
Params.ini
----------

PLEASE SEE :ref:`inifile` for detailed information about the inifile and how it is constructed

The inifile contains four subsections: filters, convergence, rand_seed, and bse. 

The `filters` subsection allows you to specify how you would like to filter the binary population. See the inifile for a description of each flag.

The `convergence` subsection allows you to specify a particular region of parameter space where you would like the convergence of each binary parameter distribution to be tested. The only implemented convergence filter is for LISA binaries, where the convergence is computed only for binaries with orbital period less than 5000 s. This allows for low probability, but high signal to noise binaries with very short orbital periods to be fully accounted for. 

The `rand_seed` subsections allows you to initialize the population with a random seed for reproduceability. Note that for each simulated binary, cosmic also returns the initial conditions, including a random seed that will reproduce that exact binary. The random seed produced in the rand_seed subsection allows full populations to be reproduced. This is particularly useful when comparing two popuations, e.g. binary black holes and binary neutron stars, which should be simulated separately, but using the same rand_seed value.

The `bse` subsection is where all the bse flags are specified.

-------------------
Sample command line
------------------- 

Below is an example command line call for runFixedPop:

.. code-block:: bash
   
   runFixedPop --final_kstar1 11 --final_kstar2 10 11 --inifile Params.ini --galaxy_component ThinDisk --metallicity 0.02 --convergence-params mass_1 mass_2 sep ecc --initial_samp multidim --Nstep 15000 --Niter 1000000000 -n 2 

A breakdown of each argument follows:

* ---final_ktar1 11 ---final_kstar2 10 11 : this tells cosmic to keep all systems where the primary star is a CO WD and the secondary star is either a CO or He WD. 

* --inifile Params.ini : this tells cosmic where to look for the BSE flags that set the binary evolution model; in this case, the inifile is assumed to be in the same directory that the command line call is made

* --galaxy_component ThinDisk : this tells cosmic to run a Milky Way thin disk, which will implement constant star formation history over 10 Gyr

* --metallicity 0.02 : this tells cosmic to use solar metallicity for every binary evolved 

* --convergence-params mass_1 mass_2 sep ecc : this tells cosmic to track how the distributions of the masses, semimajor axis, and eccentricity of the binaries change with each iterated population

* --initial_samp multidim : this tells cosmic to initialize the binaries with multidimensional parameter distributions according to `Moe & Di Stefano 2017 <http://adsabs.harvard.edu/abs/2017ApJS..230...15M>`_

* --Nstep 5000 --Niter 10000000 -n 1 : this tells cosmic to use 1 processor to evolve a maximum of 1e7 systems and check in every 5000 systems to track how the shape of the distributions of the parameters specified in convergence-params change

===================
Stopping conditions
===================

There are two stopping conditions for runFixedPop:

1. the shape of the parameter distributions for the parameters specified in --convergence-params converge to a shape, regardelss of adding new binaries to the evolved population. This is quantified by the match criteria detailed in Breivik+2018 in prep

2. the number of binaries evolved exceeds Niter. 

=====================
Output of runFixedPop
=====================

PLEASE SEE :ref:`output_info` for more information about the output data frames including
what each column means and the units.

The output of runFixed pop is the `fixed population`, an hdf5 file with a naming scheme that tells you the Galactic component and final kstars of the population; the data file created by the runFixedPop call above is: dat_ThinDisk_11_10_11.h5. 

The fixed population contains three pandas DataFrames accessed by the following keys:

* bcm : The final state of the converged population at the present epoch

* bpp : The evolutionary history of the systems in the bcm data set

* initCond : The initial conditions for each binary in the bcm data set

Each of these DataFrames shares a common 'binary_number' column which is used to index the population.


*************
gxRealization
*************
The gxRealization exectuable uses the fixed population and a model for the spatial distribution of systems in a given Galactic component to Monte Carlo sample synthetic Milky Way population realizations. The necessary arguments for the gxRealization executable can be accessed using the help:

.. code-block:: bash

   gxRealization -h 
.. program-output:: gxRealization -h 

======
Inputs
======
cosmic has several different models to spatially distribute binary sources in the Galaxy, depending on the --galaxy_component selection. These choices are detailed below, however, in `all cases` the orbital inclination, longitude of the ascending node, and the argument of periapse are randomized.

--------
ThinDisk
--------
There are three models to choose from with a ThinDisk population, where the differences between each model lie in the distbrituion of binaries above and below the disk and the scaling factor of each distribution. The radial distribution of binaries is always an exponential decay, though the scaling can vary from model to model. The azimuthal distribution is always uniform. 

* 'sech_squared' : Radial exponential decay distribution with scale factor of 2.5 kpc and sech_squared distribution with scale factor of 0.3 kpc; consistent with `Nelemans+2001 <http://adsabs.harvard.edu/abs/2001A%26A...375..890N>`_

* 'double_exp' : Radial and vertical exponential decay distributions, with a scale factor of 2.5 kpc radially and 0.3 kpc vertically

* 'McMillan' : Radial and vertical exponential decay distributions, with a scale factor of 2.9 kpc radially and 0.3 kpc vertically; consistent with `McMillan 2011 <http://adsabs.harvard.edu/abs/2011MNRAS.414.2446M>`_

-----
Bulge
-----
There are two models to choose from with a Bulge population.

* 'exp_squared' : Radial exponential squared decay distribution with a scale factor of 0.5 kpc, uniform azimuthal distribution and uniform in cos polar distribution; consistent with `Nelemans+2001 <http://adsabs.harvard.edu/abs/2001A%26A...375..890N>`_

* 'McMillan' : Three dimensional distribution consistent with `McMillan 2011 <http://adsabs.harvard.edu/abs/2011MNRAS.414.2446M>`_


---------
ThickDisk
---------
There are two models to choose from with a ThickDisk population. Both use exponentional decay distributions for the radial and vertical directions and uniform azimuthal distribution, but differ in the choice of scale factor.

* 'double_exp' : Radial and vertical exponential decay distributions with radial scale factor of 2.5 kpc and vertical scale factor of 1 kpc

* 'McMillan' : Radial and vertical exponential decay distributions with radial scale factor of 3.1 kpc and vertical scale factor of 0.9 kpc, consistent with `McMillan 2011 <http://adsabs.harvard.edu/abs/2011MNRAS.414.2446M>`_


===================
Sample command line
===================
Below is a sample command line input to run 100 Galactic realizations for a thin disk population of white dwarf binaries with a CO WD primary and CO or He WD secondary.

.. code-block:: bash
   
   gxRealization --final_kstar1 11 --final_kstar2 10 11 --galaxy_component ThinDisk --dist_model McMillan --N_realizations 100 --gx_save True --HG_save False --LISA_calc True -n 1

Let's break down each argument:

* --final_ktar1 11 --final_kstar2 10 11 : this tells cosmic that we want to keep all systems where the primary star is a CO WD and the secondary star is either a CO or He WD.

* --galaxy_component ThinDisk --dist_model McMillan : this tells cosmic to distribute the thin disk sources according to McMillan 2011

* --N_realizations 100 : this tells cosmic to generate 100 thin disk realizations

* --gx_save True : this tells cosmic to save the galactic realizations - NOTE: this can generate large amounts of data (~3.5G per realization) for large populations (e.g. white dwarf binaries) and large numbers of galactic realizations

* --HG_save False : this tells cosmic to ignore any systems that undergo a common envelope while the secondary is on the Hertzsprung Gap; these systems are expected to merge - see Belczynski et al. 2008 for more details

* --LISA_calc True : this tells cosmic to compute the gravitational wave power which can be used to generate the population's power spectral density observable by LISA

* -n : this tells cosmic to use 1 processor - NOTE: multiprocessing is available!

=======================
Output of gxRealization
=======================
The output of gxRealization is file for every Galactic realization created where each realization file is an hdf5 file with a naming scheme that tells you the Galctic realization number, the Galactic component, and final kstars of the population; the 0th realization data file created by the gxRealization call above is: gxReal_0_ThinDisk_11_10_11.h5.

Each realization contains up to four pandas DataFrames accessed by the following keys:

* gx_dat : The full realization, including binary parameters, spatial distribution, and binary orientation

* PSD : The LISA signal to power spectral density data including the gravitational wave frequency and PSD

WARNING - if you set gx_save to False and LISA_calc to False you will get no output!
