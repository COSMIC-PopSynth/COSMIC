.. _runpop:

####################################
Generate a binary population by hand
####################################

The process to generate a synthetic binary population, is similar to the
process to evolve a single/multiple binaries by hand: first generate an 
initial population, then evolve it with the Evolve class. 

An initialized binary population consists of a collection of binary systems
with assigned primary and secondary masses, orbital periods, eccentricities,
metallicities, and star formation histories. These parameters are randomly
sampled from observationally motivated distribution functions. 

In COSMIC, the initial sample is done through an initial binary sampler which works
with the InitialBinaryTable class. There are two samplers available: 

1. `independent` : initialize binaries with independent parameter 
distributions for the primary mass, mass ratio, eccentricity, separation, 
and binary fraction

2. `multidim` : initialize binaries with multidimensional parameter 
distributions according to `Moe & Di Stefano 2017 <http://adsabs.harvard.edu/abs/2017ApJS..230...15M>`_

We consider both cases below. 

***********
independent 
***********

First import the InitialBinaryTable class and the independent sampler

.. ipython::

    In [1]: from cosmic.sample.initialbinarytable import InitialBinaryTable

    In [2]: from cosmic.sample.sampler import independent

The independent sampler contains multiple models for each binary parameter.
You can access the available models using the independent sampler help call:    

.. ipython::

    In [3]: help(independent.get_independent_sampler)

The final_kstar1 and final_kstar2 parameters are lists that contain the kstar types
that you would like the final population to contain. 

The final kstar is the final state of the binary system we are interested in and is based on the BSE kstar naming convention, see :ref:`kstar-table` for more information.

Thus, if you want to generate a 
population containing double white dwarfs with CO and ONe WD primaries and He-WD secondaries, 
the final kstar inputs would be:

.. ipython::

    In [4]: final_kstar1 = [11, 12]

    In [5]: final_kstar2 = [10]

Similar to the help for the sampler, the different models that can be used for each parameter 
to be sampled can be accessed by the help function for the argument. The syntax for each parameter
sample is always: sample_`parameter`. See the example for the star formation
history (SFH) below:

.. ipython::

    In [6]: help(independent.Sample.sample_SFH)

Using the final kstar inputs above, the initial binary population is sampled as:

.. ipython::

    In [6]: InitialBinaries, mass_singles, mass_binaries, n_singles, n_binaries = InitialBinaryTable.sampler('independent', final_kstar1, final_kstar2, binfrac_model=0.5, primary_model='kroupa93', ecc_model='thermal', SFH_model='const', component_age=10000.0, met=0.02, size=10000)

    In [7]: print(InitialBinaries)

NOTE: the length of the initial binary data set, InitialBinaries, does not always match 
the size parameter provided to InitialBinaryTable.sampler. 
This is becuase the sampler accounts for a binary fraction specified by the user with the binfrac_model parameter, which is either a fraction between 0 and 1 or mass dependend following the prescription in `van Haaften+2013 <http://adsabs.harvard.edu/abs/2012A%26A...537A.104V>`_. 
Since we are interested in binaries, we only retain the binary systems that are likely to produce the user specified final kstar types. However, we also keep track of the total mass of the single and binary stars as well as the numbre of binary and single stars so that we can scale our results to larger populations. If you don't want to filter the binaries, you can supply final kstars as 

.. ipython::

    In [8]: final_kstar = range(0,15)

********
multidim
********

COSMIC implements multidimensionally distributed initial binaries according to `Moe & Di Stefano 2017 <http://adsabs.harvard.edu/abs/2017ApJS..230...15M>`_. The python code used in COSMIC to create this sample was written by Mads Sorenson, and is based on the IDL codes written to accompany `Moe & Di Stefano 2017 <http://adsabs.harvard.edu/abs/2017ApJS..230...15M>`_. 

The multidimensional initial binary data is sampled in COSMIC as follows:

.. ipython::

    In [1]: from cosmic.sample.initialbinarytable import InitialBinaryTable

    In [2]: from cosmic.sample.sampler import multidim

To see the arguments necessary to call the multidimensional sampler use the help function:

.. ipython::
 
    In [3]: help(multidim.get_multidim_sampler)  

The random seed is used to reproduce your initial sample, since there are several stochastic processes involved in the muldimensional sample. 
As in the independent sampler, the final_kstar1 and final_kstar2 inputs are lists containing the kstar types that the evolved population should contain.

The multidimensional sample is generated as follows:

.. ipython::

    In [4]: InitialBinaries, mass_singles, mass_binaries, n_singles, n_binaries = InitialBinaryTable.sampler('multidim', final_kstar1=[11], final_kstar2=[11], rand_seed=2, nproc=1, SFH_model='const', component_age=10000.0, met=0.02, size=10)

    In [5]: print(InitialBinaries)

.. note::

    NOTE that in the multidimensional case, the binary fraction is a parameter in the sample. This results in the size of the initial binary data matching the size provided to the sampler. As in the independent sampling case, we keep track of the total sampled mass of singles and binaries as well as the total number of single and binary stars to scale thesimulated population to astrophysical populations.

*************************************
Evolving an initial binary population
*************************************
As in :ref:`examples`, once an initial binary population is sampled, it is evolved using the Evolve class. Note that the same process used in :ref:`examples` applies here as well: the BSEDict must be supplied, but only need be resupplied if the flags in the dictionary change.

The syntax for the Evolve class is as follows:

.. ipython::

    In [1]: from cosmic.evolve import Evolve   

    In [2]: BSEDict = {'xi': 0.5, 'bhflag': 1, 'neta': 0.5, 'windflag': 3, 'wdflag': 0, 'alpha1': 1.0, 'pts1': 0.05, 'pts3': 0.02, 'pts2': 0.01, 'epsnov': 0.001, 'hewind': 1.0, 'ck': -1000, 'bwind': 0.0, 'lambdaf': 0.5, 'mxns': 3.0, 'beta': -1.0, 'tflag': 1, 'acc2': 1.5, 'nsflag': 4, 'ceflag': 0, 'eddfac': 1.0, 'ifflag': 0, 'bconst': -3000, 'sigma': 265.0, 'gamma': -2.0, 'pisn': 45.0, 'natal_kick_array' :[-100.0,-100.0,-100.0,-100.0,-100.0,-100.0], 'bhsigmafrac' : 1.0, 'polar_kick_angle' : 90, 'qcrit_array' : [0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0], 'cekickflag' : 2, 'cehestarflag' : 0, 'cemergeflag' : 0, 'ecsn' : 2.5, 'ecsn_mlow' : 1.4, 'aic' : 1, 'ussn' : 0, 'sigmadiv' :-20.0, 'qcflag' : 2, 'eddlimflag' : 0, 'fprimc_array' : [2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0], 'bhspinflag' : 0, 'bhspinmag' : 0.0, 'rejuv_fac' : 1.0, 'rejuvflag' : 0, 'rejuvflag' : 0, 'htpmb' : 1, 'ST_cr' : 1, 'ST_tide' : 1}

    In [3]: bpp, bcm, initC  = Evolve.evolve(initialbinarytable=InitialBinaries, BSEDict=BSEDict)

    In [4]: print(bcm.iloc[:10])

    In [5]: print(bpp)
