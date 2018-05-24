.. _examples:

########################################
Using COSMIC to evolve binaries with BSE
########################################


************
Introduction
************

COSMIC can simulate binaries for several different use cases. Below 
you'll find examples to run a single binary system, multiple binary
systems or a grid of binaries. See below for the process to simulate
a population of binaries consistent with a user-supplied star formation 
history for a single compact object population (e.g. BH-BH) or a range
of compact object populations (e.g. combinations of BH, NS, WD) as 
described in Breivik & Larson (2018).



*************************************
Creating and evolving a single binary
*************************************

Let's initialize and evolve a whopper of a binary that 
could have formed our friend GW150914

.. ipython::

    In [1]: from aCOSMIC.sample.initialbinarytable import InitialBinaryTable

    In [2]: from aCOSMIC.evolve import Evolve

    In [3]: single_binary = InitialBinaryTable.SingleBinary(m1=114.0, m2=45.0, porb=50.0, ecc=0.65, tphysf=13700.0, kstar1=1, kstar2=1, metallicity=0.002)

    In [4]: print(single_binary)

    In [5]: BSEDict = {'xi': 0.5, 'bhflag': 1, 'neta': 0.5, 'windflag': 3, 'wdflag': 0, 'alpha1': 1.0, 'pts1': 0.05, 'pts3': 0.02, 'pts2': 0.01, 'epsnov': 0.001, 'hewind': 1.0, 'CK': -1000, 'bwind': 0.0, 'lambdaf': 1.0, 'mxns': 3.0, 'beta': -1.0, 'tflag': 1, 'acc2': 1.5, 'nsflag': 3, 'ceflag': 0, 'eddfac': 1.0, 'merger': 0, 'ifflag': 0, 'bconst': -3000, 'sigma': 265.0, 'gamma': -2.0}

    In [6]: EvolvedBinaryBPP, EvolvedBinaryBCM = Evolve.evolve(initialbinarytable=single_binary, BSEDict=BSEDict)

    In [7]: print(EvolvedBinaryBPP)

    In [8]: print(EvolvedBinaryBCM)

***************************************
Creating and evolving multiple binaries
***************************************

Now let's initialize and evolve systems that could form GW150914 and GW170817

.. ipython::

    In [1]: from aCOSMIC.sample.initialbinarytable import InitialBinaryTable

    In [2]: from aCOSMIC.evolve import Evolve

    In [3]: binary_set = InitialBinaryTable.MultipleBinary(m1=[114.0, 11.8], m2=[45.0, 11.1], porb=[50.0,2211.0], ecc=[0.65,0.55], tphysf=[13700.0,13700.0], kstar1=[1,1], kstar2=[1,1], metallicity=[0.002,0.02])

    In [4]: print(binary_set)

    In [5]: EvolvedBinariesBPP, EvolvedBinariesBCM  = Evolve.evolve(initialbinarytable=binary_set, BSEDict=BSEDict)

    In [6]: print(EvolvedBinariesBPP)

    In [7]: print(EvolvedBinariesBCM)


****************************************
Creating and evolving a grid of binaries
****************************************

Sometimes it is helpful to run a grid of initial binaries to explore how
changing a single paramter affects the evolved binary. Here we will evolve 
the same system we ran for GW150914, but run over several initial orbital
periods

.. ipython::

    In [1]: from aCOSMIC.sample.initialbinarytable import InitialBinaryTable

    In [2]: from aCOSMIC.evolve import Evolve

    In [3]: n_grid = 10 

    In [4]: binary_grid = InitialBinaryTable.MultipleBinary(m1=np.ones(n_grid)*114.0, m2=np.ones(n_grid)*45.0, porb=np.logspace(0,4,n_grid), ecc=np.ones(n_grid)*0.65, tphysf=np.ones(n_grid)*13700.0, kstar1=np.ones(n_grid), kstar2=np.ones(n_grid), metallicity=np.ones(n_grid)*0.002)

    In [4]: print(binary_grid)

    In [5]: EvolvedBinariesBPP, EvolvedBinariesBCM  = Evolve.evolve(initialbinarytable=binary_grid, BSEDict=BSEDict)

    In [6]: print(EvolvedBinariesBPP)

    In [7]: print(EvolvedBinariesBCM)

*************************************************
Creating and evolving a sampled binary population
*************************************************

To generate a Milky Way like population, we need to generate an initial set of
binaries that is representative of a population born in the Milky Way. This means
we need to supply several arguments to our initial binary sampler. There are two 
distinct ways to generate the initial population:
1 - assume all binary paramters are independent of one another
2 - try to account for paramter dependencies

We consider both cases below. 

**********************************
Sampling Initial Binary Conditions
**********************************

.. ipython::

    In [1]: from aCOSMIC.sample.initialbinarytable import InitialBinaryTable

    In [2]: IBT, sampled_mass = InitialBinaryTable.sampler('independent', primary_min=0.08, primary_max=5.0, primary_model='kroupa93', ecc_model='thermal', SFH_model='const', component_age=10000.0, size=1000)

    In [3]: print(IBT)
