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

In cosmic, the initial sample is done through an initial binary sampler which works
with the InitialBinaryTable class. There are two samplers available: 

1. `independent` : initial binary parameters are distributed independently
2. `multidim` : accounts for parameter dependencies following `Moe & Di Stefano 2017 <http://adsabs.harvard.edu/abs/2017ApJS..230...15M>`_

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

The final kstar is the final state of the binary system we are interested in and is based on the BSE kstar naming conventions. The conventions are as follows:


=====     ==================================
kstar     evolutionary stage
=====     ==================================
0         MS, < 0.7 Msun
1         MS, > 0.7 Msun
2         Hertzsprung Gap
3         First Giant Branch
4         Core Helium Burning
5         Early Asymptotic Giant Branch
6         Thermally Pulsing AGB
7         Naked Helium Star MS
8         Naked Helium Star Hertzsprung Gap
9         Naked Helium Star Giant Branch
10        Helium White Dwarf
11        Carbon/Oxygen White Dwarf
12        Oxygen/Neon White Dwarf
13        Neutron Star
14        Black Hole
15        Massless Remnant
=====     ==================================

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

    In [4]: help(independent.Sample.sample_SFH)

Using the final kstar inputs above, the initial binary population is sampled as:

.. ipython::

    In [5]: InitialBinaries, sampled_mass, n_sampled = InitialBinaryTable.sampler('independent', final_kstar1, final_kstar2, primary_model='kroupa93', ecc_model='thermal', SFH_model='const', component_age=10000.0, met=0.02, size=10000)

    In [6]: print(InitialBinaries)

NOTE: the length of the initial binary data set, InitialBinaries, does not match 
the size parameter provided to InitialBinaryTable.sampler. 
This is becuase the sampler begins by sampling a set of stellar masses of size=size, then assigns each of the stellar masses to be either single or in a binary system following the prescription in `van Haaften+2013 <http://adsabs.harvard.edu/abs/2012A%26A...537A.104V>`_. 
Since we are interested in binaries, we only retain the binary systems. However, we also keep track of the total mass sampled so that we can scale our results to a full Milky Way population.

********
multidim
********

cosmic implements multidimensionally distributed initial binaries according to `Moe & Di Stefano 2017 <http://adsabs.harvard.edu/abs/2017ApJS..230...15M>`_. The python code used in cosmic to create this sample was written by Mads Sorenson, and is based on the IDL codes written to accompany `Moe & Di Stefano 2017 <http://adsabs.harvard.edu/abs/2017ApJS..230...15M>`_. 

The multidimensional initial binary data is sampled in cosmic as follows:

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

    In [5]: InitialBinaries, sampled_mass, n_sampled = InitialBinaryTable.sampler('multidim', [11], [11], 2, 1, 'const', 10000.0, 0.02, 10)

    In [6]: print(InitialBinaries)

NOTE that in the multidimensional case, the binary fraction is a parameter in the sample. This results in the size of the initial binary data matching the size provided to the sampler. As in the independent sampling case, we keep track of the total sampled mass to scale our simulated population to the full Milky Way.

*************************************
Evolving an initial binary population
*************************************
As in :ref:`examples`, once an initial binary population is sampled, it is evolved using the Evolve class. Note that the same process used in :ref:`examples` applies here as well: the BSEDict must be supplied, but only need be resupplied if the flags in the dictionary change.

The syntax for the Evolve class is as follows:

.. ipython::

    In [2]: from cosmic.evolve import Evolve   

    In [4]: BSEDict = {'xi': 0.5, 'bhflag': 1, 'neta': 0.5, 'windflag': 3, 'wdflag': 0, 'alpha1': 1.0, 'pts1': 0.05, 'pts3': 0.02, 'pts2': 0.01, 'epsnov': 0.001, 'hewind': 1.0, 'ck': -1000, 'bwind': 0.0, 'lambdaf': 1.0, 'mxns': 3.0, 'beta': -1.0, 'tflag': 1, 'acc2': 1.5, 'nsflag': 3, 'ceflag': 0, 'eddfac': 1.0, 'ifflag': 0, 'bconst': -3000, 'sigma': 265.0, 'gamma': -2.0, 'pisn': 45.0, 'natal_kick_array' :[-100.0,-100.0,-100.0,-100.0,-100.0,-100.0], 'bhsigmafrac' : 1.0, 'polar_kick_angle' : 90, 'qcrit_array' : [0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0], 'cekickflag' : 0, 'cehestarflag' : 0, 'cemergeflag' : 0, 'ecsn' : 2.5, 'ecsn_mlow' : 1.6, 'aic' : 1, 'ussn' : 0, 'sigmadiv' :-20.0, 'qcflag' : 1}

    In [5]: bpp, bcm, initC  = Evolve.evolve(initialbinarytable=InitialBinaries, BSEDict=BSEDict)

    In [6]: print(bcm.iloc[:10])

    In [7]: print(bpp)

Note that this process doesn't try to choose the `right` number of binaries to evolve. If you are interested in generating a realistic synthetic Milky Way population, you should head over to :ref:`fixedpop`. For details on the process to generate synthetic Milky Way binary populations, see Breivik+2018 (in prep). 
