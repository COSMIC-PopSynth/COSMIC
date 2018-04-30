.. _examples:

########################################
Double White Dwarf Accretion with COSMIC
########################################


************
Introduction
************

Science

*************************************
Creating and evolving a single binary
*************************************

.. ipython::

    In [1]: from aCOSMIC.sample.initialbinarytable import InitialBinaryTable

    In [2]: from aCOSMIC.evolve import Evolve

    In [3]: single_binary = InitialBinaryTable.SingleBinary(m1=1.5, m2=1.0, porb=25571606.5405, ecc=0.05, tphysf=8763.7162809763449, kstar1=1, kstar2=1, metallicity=0.02)

    In [4]: print(single_binary)

    In [5]: BSEDict = {'xi': 0.5, 'bhflag': 1, 'neta': 0.5, 'windflag': 3, 'wdflag': 0, 'alpha1': 1.0, 'pts1': 0.05, 'pts3': 0.02, 'pts2': 0.01, 'epsnov': 0.001, 'hewind': 1.0, 'CK': -1000, 'bwind': 0.0, 'lambdaf': -1.0, 'mxns': 3.0, 'beta': -1.0, 'tflag': 1, 'acc2': 1.5, 'nsflag': 2, 'ceflag': 0, 'eddfac': 1.0, 'merger': 0, 'ifflag': 0, 'bconst': -3000, 'sigma': 265.0, 'gamma': -2.0}

    In [6]: EvolvedBinaryBPP, EvolvedBinaryBCM = Evolve.evolve(initialbinarytable=single_binary, BSEDict=BSEDict)

    In [7]: print(EvolvedBinaryBPP)

    In [8]: print(EvolvedBinaryBCM)

***********************************
Creating and evolving a binary grid
***********************************

.. ipython::

    In [1]: from aCOSMIC.sample.initialbinarytable import InitialBinaryTable

    In [2]: from aCOSMIC.evolve import Evolve

    In [3]: binary_grid = InitialBinaryTable.MultipleBinary(m1=[1.5, 2.0], m2=[1.0, 1.75], porb=[2557160646789.5405,2557160646789.5405], ecc=[0.05,0.01], tphysf=[8763.7162809763449,10000.7162809763449], kstar1=[11,12], kstar2=[11,11], metallicity=[0.02,0.02])

    In [4]: print(binary_grid)

    In [5]: EvolvedBinariesBPP, EvolvedBinariesBCM  = Evolve.evolve(initialbinarytable=binary_grid, BSEDict=BSEDict)

    In [6]: print(EvolvedBinariesBPP)

    In [7]: print(EvolvedBinariesBCM)

*************************************************
Creating and evolving a sampled binary population
*************************************************

.. ipython::

    In [1]: from aCOSMIC.sample.initialbinarytable import InitialBinaryTable

    In [2]: sampled_binaries = InitialBinaryTable.sampler('independent')

    In [3]: print(sampled_binaries)

**********************************
Sampling Initial Binary Conditions
**********************************
