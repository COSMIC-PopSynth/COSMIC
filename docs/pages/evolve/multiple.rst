**************************
Evolving multiple binaries
**************************

Following on from `evolving a single binary <single.html>`_, we can also evolve multiple binaries.
Let's start by importing the necessary modules:

.. ipython::

    In [1]: from cosmic.sample.initialbinarytable import InitialBinaryTable

    In [2]: from cosmic.evolve import Evolve


And use the same BSE dict as before:

.. ipython::

    In [3]: BSEDict = {'xi': 1.0, 'bhflag': 1, 'neta': 0.5, 'windflag': 3, 'wdflag': 1, 'alpha1': 1.0, 'pts1': 0.001, 'pts3': 0.02, 'pts2': 0.01, 'epsnov': 0.001, 'hewind': 0.5, 'ck': 1000, 'bwind': 0.0, 'lambdaf': 0.0, 'mxns': 3.0, 'beta': -1.0, 'tflag': 1, 'acc2': 1.5, 'grflag' : 1, 'remnantflag': 4, 'ceflag': 0, 'eddfac': 1.0, 'ifflag': 0, 'bconst': 3000, 'sigma': 265.0, 'gamma': -2.0, 'pisn': 45.0, 'natal_kick_array' : [[-100.0,-100.0,-100.0,-100.0,0.0], [-100.0,-100.0,-100.0,-100.0,0.0]], 'bhsigmafrac' : 1.0, 'polar_kick_angle' : 90, 'qcrit_array' : [0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0], 'cekickflag' : 2, 'cehestarflag' : 0, 'cemergeflag' : 0, 'ecsn' : 2.25, 'ecsn_mlow' : 1.6, 'aic' : 1, 'ussn' : 0, 'sigmadiv' :-20.0, 'qcflag' : 1, 'eddlimflag' : 0, 'fprimc_array' : [2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0], 'bhspinflag' : 0, 'bhspinmag' : 0.0, 'rejuv_fac' : 1.0, 'rejuvflag' : 0, 'htpmb' : 1, 'ST_cr' : 1, 'ST_tide' : 1, 'bdecayfac' : 1, 'rembar_massloss' : 0.5, 'kickflag' : 1, 'zsun' : 0.014, 'bhms_coll_flag' : 0, 'don_lim' : -1, 'acc_lim' : -1, 'rtmsflag' : 0, 'wd_mass_lim': 1}


Below is an example for systems that could form GW150914 and GW170817 - like binaries.

.. ipython::
    :okwarning:

    In [3]: binary_set = InitialBinaryTable.InitialBinaries(m1=[85.543645, 11.171469], m2=[84.99784, 6.67305],
       ...:                                                 porb=[446.795757, 170.758343], ecc=[0.448872, 0.370],
       ...:                                                 tphysf=[13700.0, 13700.0],
       ...:                                                 kstar1=[1, 1], kstar2=[1, 1],
       ...:                                                 metallicity=[0.002, 0.02])

    In [4]: print(binary_set)

    In [5]: import numpy as np

    In [6]: np.random.seed(5)

    In [7]: bpp, bcm, initC, kick_info = Evolve.evolve(initialbinarytable=binary_set, BSEDict=BSEDict)

As before, bpp, bcm, and initC are returned as pandas DataFrames which assign an
index to each binary system we evolve. We can access each binary as follows:

.. ipython::

    In [8]: print(bpp.loc[0])

    In [9]: print(bcm.loc[0])

    In [10]: print(initC.loc[0])

    In [11]: print(bpp.loc[1])

The plotting function can also take in multiple binaries. Let's plot both the GW150914-like
progenitor evolution and the GW170817-like progenitor evolutions. For the GW170817-like
progenitor, we expect most of the evolution to take place in the first ~60 Myr.

.. ipython::
    :okwarning:
    :okexcept:

    In [12]: fig = evolve_and_plot(binary_set, t_min=None, t_max=[6.0, 60.0], BSEDict=BSEDict, sys_obs={})


.. plot::

    from cosmic.sample.initialbinarytable import InitialBinaryTable
    from cosmic.plotting import evolve_and_plot
    import numpy as np
    np.random.seed(5)
    BSEDict = {'xi': 1.0, 'bhflag': 1, 'neta': 0.5, 'windflag': 3, 'wdflag': 1, 'alpha1': 1.0, 'pts1': 0.001, 'pts3': 0.02, 'pts2': 0.01, 'epsnov': 0.001, 'hewind': 0.5, 'ck': 1000, 'bwind': 0.0, 'lambdaf': 0.0, 'mxns': 3.0, 'beta': -1.0, 'tflag': 1, 'acc2': 1.5, 'grflag' : 1, 'remnantflag': 4, 'ceflag': 0, 'eddfac': 1.0, 'ifflag': 0, 'bconst': 3000, 'sigma': 265.0, 'gamma': -2.0, 'pisn': 45.0, 'natal_kick_array' : [[-100.0,-100.0,-100.0,-100.0,0.0], [-100.0,-100.0,-100.0,-100.0,0.0]], 'bhsigmafrac' : 1.0, 'polar_kick_angle' : 90, 'qcrit_array' : [0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0], 'cekickflag' : 2, 'cehestarflag' : 0, 'cemergeflag' : 0, 'ecsn' : 2.25, 'ecsn_mlow' : 1.6, 'aic' : 1, 'ussn' : 0, 'sigmadiv' :-20.0, 'qcflag' : 1, 'eddlimflag' : 0, 'fprimc_array' : [2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0], 'bhspinflag' : 0, 'bhspinmag' : 0.0, 'rejuv_fac' : 1.0, 'rejuvflag' : 0, 'htpmb' : 1, 'ST_cr' : 1, 'ST_tide' : 1, 'bdecayfac' : 1, 'rembar_massloss' : 0.5, 'kickflag' : 1, 'zsun' : 0.014, 'bhms_coll_flag' : 0, 'don_lim' : -1, 'acc_lim' : -1, 'rtmsflag' : 0, 'wd_mass_lim': 1}
    binary_set = InitialBinaryTable.InitialBinaries(m1=[85.543645, 11.171469], m2=[84.99784, 9.67305], porb=[446.795757, 370.758343], ecc=[0.448872, 0.370], tphysf=[13700.0, 13700.0], kstar1=[1, 1], kstar2=[1, 1], metallicity=[0.002, 0.02])
    fig = evolve_and_plot(binary_set, t_min=None, t_max=[6.0, 60.0], BSEDict=BSEDict, sys_obs={})

