.. _runpop:

##################################
Evolve a binary population by hand
##################################

To create a synthetic binary population, we need to generate an initial set 
of binaries that is representative of a population born in the Milky Way. 
This means we need to supply several arguments to our initial binary sampler. 
There are two  ways to generate the initial population:

1. assume binary parameters are independent of one another
2. account for parameter dependencies following Moe & Di Stefano 2017

We consider both cases below. 

************************************************************************
Sampling an initial population with independently distributed parameters
************************************************************************

There are several different initial parameter distributions that are used 
in binary population synthesis codes. COSMIC is equipped with different 
models for each binary parameter. You can access the available models using
the help for the independent sampler as shown in the example below 

.. ipython::

    In [1]: from cosmic.sample.initialbinarytable import InitialBinaryTable

    In [2]: from cosmic.sample.sampler import independent

To see the arguments the initial binary sampler requires, use the help function

.. ipython::

    In [3]: help(independent.get_independent_sampler)

Similarly, to see the different models that can be used for each parameter 
sample, use the help function for the argument. The syntax for each parameter
sample is always: sample_`paramter`. See the example for the star formation
history (SFH) below:

.. ipython::

    In [4]: help(independent.Sample.sample_SFH)

Now let's create our first initial binary sample:

.. ipython::

    In [5]: InitialBinaries, sampled_mass = InitialBinaryTable.sampler('independent', 0.08, 5.0, 'kroupa93', 'thermal', 'const', 10000.0, 0.02, 10)

    In [6]: print(InitialBinaries)

NOTE: the length of the initial binary data set, IBT, does not match 
the size parameter provided to InitialBinaryTable.sampler. 
This is becuase the sampler begins by sampling a set of stellar masses of size=size, then assigns each of the stellar masses to be either single or in a binary system. Since we are interested in binaries, we only retain the binary systems. However, we also keep track of the total mass sampled so that we can scale our findings to a full Milky Way population.

****************************************************************************
Sampling an initial population with multidimensional parameter distributions
****************************************************************************

cosmic implements multidimensionally distributed initial binaries according to Moe & Di Stefano 2017. The python code used in cosmic to create this sample was lightly adapted from python code written by Mads Sorenson, which is based on the IDL scripts written to accompany Moe & Di Stefano 2017. 

The multidimensional initial binary data is sampled in cosmic as follows:

.. ipython::

    In [1]: from cosmic.sample.initialbinarytable import InitialBinaryTable

    In [2]: from cosmic.sample.sampler import multidim

To see the arguments necessary to call the multidimensional sampler use the help function:

.. ipython::
 
    In [3]: help(multidim.get_multidim_sampler)  

The random seed is used to reproduce your initial sample, since there are several stochastic processes involved in the muldimensional sample. 

Now let's create a multidimensional initial binary sample:

.. ipython::

    In [5]: InitialBinaries, sampled_mass = InitialBinaryTable.sampler('multidim', ['11'], ['11'], 2, 1, 'const', 10000.0, 0.02, 10)

    In [6]: print(InitialBinaries)

NOTE that in the multidimensional case, the binary fraction is one of the dependent parameters. This results in the size of the initial binary data matching the size provided to the sampler. As in the independent sampling case, we keep track of the total sampled mass to scale our simulated population to the full Milky Way.

Also not that instead of supplying a minimum or maximum primary mass, we specified the final kstars. The final kstar is the final state of the binary system we are interested in and is based on the BSE kstar naming conventions. The conventions are as follows:

*   0 :        MS, < 0.7 Msun
*   1 :        MS, > 0.7 Msun
*   2 :        Hertzsprung Gap
*   3 :        First Giant Branch
*   4 :        Core Helium Burning
*   5 :        Early Asymptotic Giant Branch
*   6 :        Thermally Pulsing AGB
*   7 :        Naked Helium Star MS
*   8 :        Naked Helium Star Hertzsprung Gap
*   9 :        Naked Helium Star Giant Branch
*  10 :        Helium White Dwarf
*  11 :        Carbon/Oxygen White Dwarf
*  12 :        Oxygen/Neon White Dwarf
*  13 :        Neutron Star
*  14 :        Black Hole
*  15 :        Massless Remnant

***********************************************************
Evolving an initial binary population with the Evolve class
***********************************************************
As in :ref:`examples`, now that we have an initial binary population, we can simply evolve it using the Evolve class. The syntax is as follows:

.. ipython::

    In [2]: from cosmic.evolve import Evolve   

    In [4]: BSEDict = {'xi': 0.5, 'bhflag': 1, 'neta': 0.5, 'windflag': 3, 'wdflag': 0, 'alpha1': 1.0, 'pts1': 0.05, 'pts3': 0.02, 'pts2': 0.01, 'epsnov': 0.001, 'hewind': 1.0, 'CK': -1000, 'bwind': 0.0, 'lambdaf': 1.0, 'mxns': 3.0, 'beta': -1.0, 'tflag': 1, 'acc2': 1.5, 'nsflag': 3, 'ceflag': 0, 'eddfac': 1.0, 'merger': 0, 'ifflag': 0, 'bconst': -3000, 'sigma': 265.0, 'gamma': -2.0}

    In [5]: EvolvedBinariesBPP, EvolvedBinariesBCM, initialConditions  = Evolve.evolve(initialbinarytable=InitialBinaries, BSEDict=BSEDict)

    In [6]: print(EvolvedBinariesBCM.iloc[:10])

    In [7]: print(EvolvedBinariesBPP)

The BPP and the BCM arrays are named to follow the BSE convention. The EvolvedBinariesBPP DataFrame contains the evolutionary history of the binary and it's paramters. The EvolvedBinariesBCM DataFrame contains the current state of thebinaries at the present epoch.

Note that this process doesn't try to choose the `right` number of binaries to evolve. If you are interested in generating a realistic synthetic Milky Way population, you should head over to :ref:`fixedpop`. For details on the process to generate synthetic Milky Way binary populations, see Breivik et al 2018 (in prep). 
