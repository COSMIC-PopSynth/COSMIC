.. _examples:

###############################
Using COSMIC to evolve binaries
###############################


COSMIC can evolve binaries for several different use cases. Below
you'll find examples to run a single binary system, multiple binary
systems or a grid of binaries.


*************
single binary
*************

Below is the process to initialize and evolve a binary that
could have formed a GW150914-like binary. First, import the modules in COSMIC
that initialize and evolve the binary.

.. ipython::

    In [1]: from cosmic.sample.initialbinarytable import InitialBinaryTable

    In [2]: from cosmic.evolve import Evolve


To initialize a single binary, populate the InitialBinaries method in the
InitialBinaryTable class. Each initialized binary requires the following parameters:


* m1 : ZAMS mass of the primary star in :math:`M_{\odot}`

* m2 : ZAMS mass of the secondary star in :math:`M_{\odot}`

* porb : initial orbital period in days

* ecc : initial eccentricity

* tphysf : total evolution time of the binary in Myr

* kstar1 : initial primary stellar type, following the BSE convention

* kstar2 : initial secondary stellar type, following the BSE convention

* metallicity : metallicity of the population (e.g. :math:`Z_{\odot}=0.014`)

.. ipython::

    In [3]: single_binary = InitialBinaryTable.InitialBinaries(m1=85.543645, m2=84.99784, porb=446.795757, ecc=0.448872, tphysf=13700.0, kstar1=1, kstar2=1, metallicity=0.002)

    In [4]: print(single_binary)


The flags for the various binary evolution prescriptions used in BSE also need to be set.
Each flag is saved in the BSEDict dictionary. Note that the BSEDict
only needs to be specified the first time a binary is evolved with COSMIC or
if you need to change the binary evolution prescriptions.

If you are unfamiliar with these prescriptions, it is highly
advised to either run the defaults from the COSMIC install which are consistent
with `Breivik+2020 <https://ui.adsabs.harvard.edu/abs/2019arXiv191100903B/abstract>`_

.. ipython::

    In [5]: BSEDict = {'xi': 1.0, 'bhflag': 1, 'neta': 0.5, 'windflag': 3, 'wdflag': 1, 'alpha1': 1.0, 'pts1': 0.001, 'pts3': 0.02, 'pts2': 0.01, 'epsnov': 0.001, 'hewind': 0.5, 'ck': 1000, 'bwind': 0.0, 'lambdaf': 0.0, 'mxns': 3.0, 'beta': -1.0, 'tflag': 1, 'acc2': 1.5, 'grflag' : 1, 'remnantflag': 4, 'ceflag': 0, 'eddfac': 1.0, 'ifflag': 0, 'bconst': 3000, 'sigma': 265.0, 'gamma': -2.0, 'pisn': 45.0, 'natal_kick_array' : [[-100.0,-100.0,-100.0,-100.0,0.0], [-100.0,-100.0,-100.0,-100.0,0.0]], 'bhsigmafrac' : 1.0, 'polar_kick_angle' : 90, 'qcrit_array' : [0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0], 'cekickflag' : 2, 'cehestarflag' : 0, 'cemergeflag' : 0, 'ecsn' : 2.25, 'ecsn_mlow' : 1.6, 'aic' : 1, 'ussn' : 0, 'sigmadiv' :-20.0, 'qcflag' : 1, 'eddlimflag' : 0, 'fprimc_array' : [2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0], 'bhspinflag' : 0, 'bhspinmag' : 0.0, 'rejuv_fac' : 1.0, 'rejuvflag' : 0, 'htpmb' : 1, 'ST_cr' : 1, 'ST_tide' : 1, 'bdecayfac' : 1, 'rembar_massloss' : 0.5, 'kickflag' : 0, 'zsun' : 0.014, 'bhms_coll_flag' : 0, 'don_lim' : -1, 'acc_lim' : -1, 'rtmsflag' : 0}


Once the binary is initialized and the BSE model is set, the system is evolved with the
the Evolve class, which calls the evolv2.f subroutine in the BSE source code.

.. ipython::
    :okwarning:

    In [6]: bpp, bcm, initC, kick_info = Evolve.evolve(initialbinarytable=single_binary, BSEDict=BSEDict)


For every evolved binary system, BSE generates two arrays, which are stored as pandas DataFrames in COSMIC:

* bpp - contains binary parameters at important stages in the binary's evolution, including stellar evolutionary phase changes or mass transfer episodes.

* bcm - contains several binary parameters at user specified time steps during the binary's evolution. The default setting in COSMIC is to output the final stage of the binary at the evolution time specified by the user.

You can see the different parameters included in each DataFrame using the columns attribute of the DataFrame:

.. ipython::

    In [7]: print(bpp.columns)

    In [8]: print(bcm.columns)


The units are broadly consistent with BSE and are described in :ref:`output_info`.

The evol_type column in bpp indicates the evolutionary change that occurred for each line.
The meaning of each number is described here, :ref:`evolve-type-table`.

Each of the parameters in bpp or bcm can be accessed in the usual way for DataFrames:

.. ipython::

    In [9]: bpp.mass_1

    In [10]: bpp = bpp[['mass_1', 'mass_2', 'kstar_1', 'kstar_2', 'sep', 'evol_type']]


You can use the ``utils.convert_kstar_evol_type`` function to convert the
``kstar_1``, ``kstar_2``, and ``evol_type`` columns from integers to strings
that describe each int:

.. ipython::

    In [11]: from cosmic.utils import convert_kstar_evol_type

    In [12]: convert_kstar_evol_type(bpp)


Note that ``utils.convert_kstar_evol_type`` is only applicable to the bpp
array.

You can also use the built in plotting function to see how the system evolves:

.. ipython::
    :okwarning:

    In [12]: from cosmic.plotting import evolve_and_plot

    In [13]: single_binary = InitialBinaryTable.InitialBinaries(m1=85.543645, m2=84.99784, porb=446.795757, ecc=0.448872, tphysf=13700.0, kstar1=1, kstar2=1, metallicity=0.002)

    In [14]: BSEDict = {'xi': 1.0, 'bhflag': 1, 'neta': 0.5, 'windflag': 3, 'wdflag': 1, 'alpha1': 1.0, 'pts1': 0.001, 'pts3': 0.02, 'pts2': 0.01, 'epsnov': 0.001, 'hewind': 0.5, 'ck': 1000, 'bwind': 0.0, 'lambdaf': 0.0, 'mxns': 3.0, 'beta': -1.0, 'tflag': 1, 'acc2': 1.5, 'grflag' : 1, 'remnantflag': 4, 'ceflag': 0, 'eddfac': 1.0, 'ifflag': 0, 'bconst': 3000, 'sigma': 265.0, 'gamma': -2.0, 'pisn': 45.0, 'natal_kick_array' : [[-100.0,-100.0,-100.0,-100.0,0.0], [-100.0,-100.0,-100.0,-100.0,0.0]], 'bhsigmafrac' : 1.0, 'polar_kick_angle' : 90, 'qcrit_array' : [0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0], 'cekickflag' : 2, 'cehestarflag' : 0, 'cemergeflag' : 0, 'ecsn' : 2.25, 'ecsn_mlow' : 1.6, 'aic' : 1, 'ussn' : 0, 'sigmadiv' :-20.0, 'qcflag' : 1, 'eddlimflag' : 0, 'fprimc_array' : [2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0], 'bhspinflag' : 0, 'bhspinmag' : 0.0, 'rejuv_fac' : 1.0, 'rejuvflag' : 0, 'htpmb' : 1, 'ST_cr' : 1, 'ST_tide' : 1, 'bdecayfac' : 1, 'rembar_massloss' : 0.5, 'kickflag' : 0, 'zsun' : 0.014, 'bhms_coll_flag' : 0, 'don_lim' : -1, 'acc_lim' : -1, 'rtmsflag' : 0}

    In [15]: fig = evolve_and_plot(single_binary, t_min=None, t_max=None, BSEDict=BSEDict, sys_obs={})

.. plot::

    from cosmic.sample.initialbinarytable import InitialBinaryTable
    from cosmic.plotting import evolve_and_plot
    single_binary = InitialBinaryTable.InitialBinaries(m1=85.543645, m2=84.99784, porb=446.795757, ecc=0.448872, tphysf=13700.0, kstar1=1, kstar2=1, metallicity=0.002)
    BSEDict = {'xi': 1.0, 'bhflag': 1, 'neta': 0.5, 'windflag': 3, 'wdflag': 1, 'alpha1': 1.0, 'pts1': 0.001, 'pts3': 0.02, 'pts2': 0.01, 'epsnov': 0.001, 'hewind': 0.5, 'ck': 1000, 'bwind': 0.0, 'lambdaf': 0.0, 'mxns': 3.0, 'beta': -1.0, 'tflag': 1, 'acc2': 1.5, 'grflag' : 1, 'remnantflag': 4, 'ceflag': 0, 'eddfac': 1.0, 'ifflag': 0, 'bconst': 3000, 'sigma': 265.0, 'gamma': -2.0, 'pisn': 45.0, 'natal_kick_array' : [[-100.0,-100.0,-100.0,-100.0,0.0], [-100.0,-100.0,-100.0,-100.0,0.0]], 'bhsigmafrac' : 1.0, 'polar_kick_angle' : 90, 'qcrit_array' : [0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0], 'cekickflag' : 2, 'cehestarflag' : 0, 'cemergeflag' : 0, 'ecsn' : 2.25, 'ecsn_mlow' : 1.6, 'aic' : 1, 'ussn' : 0, 'sigmadiv' :-20.0, 'qcflag' : 1, 'eddlimflag' : 0, 'fprimc_array' : [2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0], 'bhspinflag' : 0, 'bhspinmag' : 0.0, 'rejuv_fac' : 1.0, 'rejuvflag' : 0, 'htpmb' : 1, 'ST_cr' : 1, 'ST_tide' : 1, 'bdecayfac' : 1, 'rembar_massloss' : 0.5, 'kickflag' : 0, 'zsun' : 0.014, 'bhms_coll_flag' : 0, 'don_lim' : -1, 'acc_lim' : -1, 'rtmsflag' : 0}
    fig = evolve_and_plot(single_binary, t_min=None, t_max=None, BSEDict=BSEDict, sys_obs={})


In this case, all the action happens in the first few Myr, so let's specify a t_max:

.. ipython::
    :okwarning:

    In [13]: fig = evolve_and_plot(initC, t_min=None, t_max=6.0, BSEDict={}, sys_obs={})

.. plot::

    from cosmic.sample.initialbinarytable import InitialBinaryTable
    from cosmic.plotting import evolve_and_plot
    single_binary = InitialBinaryTable.InitialBinaries(m1=85.543645, m2=84.99784, porb=446.795757, ecc=0.448872, tphysf=13700.0, kstar1=1, kstar2=1, metallicity=0.002)
    BSEDict = {'xi': 1.0, 'bhflag': 1, 'neta': 0.5, 'windflag': 3, 'wdflag': 1, 'alpha1': 1.0, 'pts1': 0.001, 'pts3': 0.02, 'pts2': 0.01, 'epsnov': 0.001, 'hewind': 0.5, 'ck': 1000, 'bwind': 0.0, 'lambdaf': 0.0, 'mxns': 3.0, 'beta': -1.0, 'tflag': 1, 'acc2': 1.5, 'grflag' : 1, 'remnantflag': 4, 'ceflag': 0, 'eddfac': 1.0, 'ifflag': 0, 'bconst': 3000, 'sigma': 265.0, 'gamma': -2.0, 'pisn': 45.0, 'natal_kick_array' : [[-100.0,-100.0,-100.0,-100.0,0.0], [-100.0,-100.0,-100.0,-100.0,0.0]], 'bhsigmafrac' : 1.0, 'polar_kick_angle' : 90, 'qcrit_array' : [0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0], 'cekickflag' : 2, 'cehestarflag' : 0, 'cemergeflag' : 0, 'ecsn' : 2.25, 'ecsn_mlow' : 1.6, 'aic' : 1, 'ussn' : 0, 'sigmadiv' :-20.0, 'qcflag' : 1, 'eddlimflag' : 0, 'fprimc_array' : [2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0], 'bhspinflag' : 0, 'bhspinmag' : 0.0, 'rejuv_fac' : 1.0, 'rejuvflag' : 0, 'htpmb' : 1, 'ST_cr' : 1, 'ST_tide' : 1, 'bdecayfac' : 1, 'rembar_massloss' : 0.5, 'kickflag' : 0, 'zsun' : 0.014, 'bhms_coll_flag' : 0, 'don_lim' : -1, 'acc_lim' : -1, 'rtmsflag' : 0}
    fig = evolve_and_plot(single_binary, t_min=None, t_max=6.0, BSEDict=BSEDict, sys_obs={})

*****************
multiple binaries
*****************

Multiple systems can also be initialized and evolved; below is an example for systems
that could form GW150914 and GW170817 - like binaries.

.. ipython::
    :okwarning:

    In [11]: binary_set = InitialBinaryTable.InitialBinaries(m1=[85.543645, 11.171469], m2=[84.99784, 6.67305], porb=[446.795757, 170.758343], ecc=[0.448872, 0.370], tphysf=[13700.0, 13700.0], kstar1=[1, 1], kstar2=[1, 1], metallicity=[0.002, 0.02])

    In [12]: print(binary_set)

    In [14]: import numpy as np

    In [15]: np.random.seed(5)

    In [13]: bpp, bcm, initC, kick_info = Evolve.evolve(initialbinarytable=binary_set, BSEDict=BSEDict)

Note that the BSEDict did not be reinitialized since the BSE model did not change.

As before, bpp, bcm, and initC are returned as pandas DataFrames which assign an
index to each binary system we evolve. We can access each binary as follows:

.. ipython::

    In [14]: print(bpp.loc[0])

    In [15]: print(bcm.loc[0])

    In [16]: print(initC.loc[0])

    In [17]: print(bpp.loc[1])

The plotting function can also take in multiple binaries. Let's plot both the GW150914-like
progenitor evolution and the GW170817-like progenitor evolutions. For the GW170817-like
progenitor, we expect most of the evolution to take place in the first ~60 Myr.

.. ipython::
    :okwarning:
    :okexcept:

    In [14]: fig = evolve_and_plot(binary_set, t_min=None, t_max=[6.0, 60.0], BSEDict=BSEDict, sys_obs={})


.. plot::

    from cosmic.sample.initialbinarytable import InitialBinaryTable
    from cosmic.plotting import evolve_and_plot
    import numpy as np
    np.random.seed(5)
    BSEDict = {'xi': 1.0, 'bhflag': 1, 'neta': 0.5, 'windflag': 3, 'wdflag': 1, 'alpha1': 1.0, 'pts1': 0.001, 'pts3': 0.02, 'pts2': 0.01, 'epsnov': 0.001, 'hewind': 0.5, 'ck': 1000, 'bwind': 0.0, 'lambdaf': 0.0, 'mxns': 3.0, 'beta': -1.0, 'tflag': 1, 'acc2': 1.5, 'grflag' : 1, 'remnantflag': 4, 'ceflag': 0, 'eddfac': 1.0, 'ifflag': 0, 'bconst': 3000, 'sigma': 265.0, 'gamma': -2.0, 'pisn': 45.0, 'natal_kick_array' : [[-100.0,-100.0,-100.0,-100.0,0.0], [-100.0,-100.0,-100.0,-100.0,0.0]], 'bhsigmafrac' : 1.0, 'polar_kick_angle' : 90, 'qcrit_array' : [0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0], 'cekickflag' : 2, 'cehestarflag' : 0, 'cemergeflag' : 0, 'ecsn' : 2.25, 'ecsn_mlow' : 1.6, 'aic' : 1, 'ussn' : 0, 'sigmadiv' :-20.0, 'qcflag' : 1, 'eddlimflag' : 0, 'fprimc_array' : [2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0], 'bhspinflag' : 0, 'bhspinmag' : 0.0, 'rejuv_fac' : 1.0, 'rejuvflag' : 0, 'htpmb' : 1, 'ST_cr' : 1, 'ST_tide' : 1, 'bdecayfac' : 1, 'rembar_massloss' : 0.5, 'kickflag' : 0, 'zsun' : 0.014, 'bhms_coll_flag' : 0, 'don_lim' : -1, 'acc_lim' : -1, 'rtmsflag' : 0}
    binary_set = InitialBinaryTable.InitialBinaries(m1=[85.543645, 11.171469], m2=[84.99784, 9.67305], porb=[446.795757, 370.758343], ecc=[0.448872, 0.370], tphysf=[13700.0, 13700.0], kstar1=[1, 1], kstar2=[1, 1], metallicity=[0.002, 0.02])
    fig = evolve_and_plot(binary_set, t_min=None, t_max=[6.0, 60.0], BSEDict=BSEDict, sys_obs={})


****************
grid of binaries
****************

Sometimes it is helpful to run a grid of initial binaries to explore how
changing a single paramter affects the evolved binary. Here we evolve
the same system that produces a GW150914-like binary, but run over several initial orbital
periods spaced evenly in log space.

.. ipython::
    :okwarning:

    In [16]: n_grid = 10

    In [17]: binary_grid = InitialBinaryTable.InitialBinaries(m1=np.ones(n_grid)*100.0, m2=np.ones(n_grid)*85.0, porb=np.logspace(3,5,n_grid), ecc=np.ones(n_grid)*0.65, tphysf=np.ones(n_grid)*13700.0, kstar1=np.ones(n_grid), kstar2=np.ones(n_grid), metallicity=np.ones(n_grid)*0.005)

    In [18]: BSEDict = {'xi': 1.0, 'bhflag': 1, 'neta': 0.5, 'windflag': 3, 'wdflag': 1, 'alpha1': 1.0, 'pts1': 0.001, 'pts3': 0.02, 'pts2': 0.01, 'epsnov': 0.001, 'hewind': 0.5, 'ck': 1000, 'bwind': 0.0, 'lambdaf': 0.0, 'mxns': 3.0, 'beta': -1.0, 'tflag': 1, 'acc2': 1.5, 'grflag' : 1, 'remnantflag': 4, 'ceflag': 0, 'eddfac': 1.0, 'ifflag': 0, 'bconst': 3000, 'sigma': 265.0, 'gamma': -2.0, 'pisn': 45.0, 'natal_kick_array' : [[-100.0,-100.0,-100.0,-100.0,0.0], [-100.0,-100.0,-100.0,-100.0,0.0]], 'bhsigmafrac' : 1.0, 'polar_kick_angle' : 90, 'qcrit_array' : [0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0], 'cekickflag' : 2, 'cehestarflag' : 0, 'cemergeflag' : 0, 'ecsn' : 2.25, 'ecsn_mlow' : 1.6, 'aic' : 1, 'ussn' : 0, 'sigmadiv' :-20.0, 'qcflag' : 1, 'eddlimflag' : 0, 'fprimc_array' : [2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0], 'bhspinflag' : 0, 'bhspinmag' : 0.0, 'rejuv_fac' : 1.0, 'rejuvflag' : 0, 'htpmb' : 1, 'ST_cr' : 1, 'ST_tide' : 1, 'bdecayfac' : 1, 'rembar_massloss' : 0.5, 'kickflag' : 0, 'zsun' : 0.014, 'bhms_coll_flag' : 0, 'don_lim' : -1, 'acc_lim' : -1, 'rtmsflag' : 0}


    In [18]: print(binary_grid)

    In [19]: bpp, bcm, initC, kick_info = Evolve.evolve(initialbinarytable=binary_grid, BSEDict=BSEDict)

    In [20]: print(bpp)

    In [21]: print(bcm)


*********************************************
dynamically set time resolution for bcm array
*********************************************

COSMIC has the ability to set time resolution of the bcm array depending on the current state of the evolution. Below we demonstrate three scenarios, setting dtp only during mass transfer, setting dtp to the same resolution for all of the evolution except for after the system merges or is disrupted, and finally an example of setting dtp only during the HMB stage of the evolution.

First, print all time steps during mass transfer

.. ipython::
    :okwarning:

    In [16]: single_binary = InitialBinaryTable.InitialBinaries(m1=7.806106, m2=5.381412, porb=2858.942021, ecc=0.601408, tphysf=13700.0, kstar1=1, kstar2=1, metallicity=0.02)

    In [16]: BSEDict = {'xi': 1.0, 'bhflag': 1, 'neta': 0.5, 'windflag': 3, 'wdflag': 1, 'alpha1': 1.0, 'pts1': 0.001, 'pts3': 0.02, 'pts2': 0.01, 'epsnov': 0.001, 'hewind': 0.5, 'ck': 1000, 'bwind': 0.0, 'lambdaf': 0.0, 'mxns': 3.0, 'beta': -1.0, 'tflag': 1, 'acc2': 1.5, 'grflag' : 1, 'remnantflag': 4, 'ceflag': 0, 'eddfac': 1.0, 'ifflag': 0, 'bconst': 3000, 'sigma': 265.0, 'gamma': -2.0, 'pisn': 45.0, 'natal_kick_array' : [[-100.0,-100.0,-100.0,-100.0,0.0], [-100.0,-100.0,-100.0,-100.0,0.0]], 'bhsigmafrac' : 1.0, 'polar_kick_angle' : 90, 'qcrit_array' : [0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0], 'cekickflag' : 2, 'cehestarflag' : 0, 'cemergeflag' : 0, 'ecsn' : 2.25, 'ecsn_mlow' : 1.6, 'aic' : 1, 'ussn' : 0, 'sigmadiv' :-20.0, 'qcflag' : 1, 'eddlimflag' : 0, 'fprimc_array' : [2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0], 'bhspinflag' : 0, 'bhspinmag' : 0.0, 'rejuv_fac' : 1.0, 'rejuvflag' : 0, 'htpmb' : 1, 'ST_cr' : 1, 'ST_tide' : 1, 'bdecayfac' : 1, 'rembar_massloss' : 0.5, 'kickflag' : 0, 'zsun' : 0.014, 'bhms_coll_flag' : 0, 'don_lim' : -1, 'acc_lim' : -1, 'rtmsflag' : 0}

    In [16]: bpp, bcm, initC, kick_info = Evolve.evolve(initialbinarytable=single_binary, BSEDict=BSEDict, timestep_conditions =[['RRLO_1>=1', 'dtp=0.0'], ['RRLO_2>=1', 'dtp=0.0']])

    In [16]: print(bcm[['tphys', 'kstar_1', 'kstar_2', 'mass_1', 'mass_2', 'RRLO_1', 'RRLO_2']])

Second, pick a certain resolution for the bcm array until the system mergers or is disrutped and then only print the final state

.. ipython::
    :okwarning:

    In [16]: bpp, bcm, initC, kick_info = Evolve.evolve(initialbinarytable=single_binary, BSEDict=BSEDict, timestep_conditions =[['binstate=0', 'dtp=1.0']])

    In [16]: print(bcm[['tphys', 'kstar_1', 'kstar_2', 'mass_1', 'mass_2', 'bin_state']])

Finally, we show how to print a fine resolution only during the HMXB stage of the evolution.

.. ipython::
    :okwarning:

    In [3]: single_binary = InitialBinaryTable.InitialBinaries(m1=85.543645, m2=84.99784, porb=446.795757, ecc=0.448872, tphysf=13700.0, kstar1=1, kstar2=1, metallicity=0.002)

    In [5]: BSEDict = {'xi': 1.0, 'bhflag': 1, 'neta': 0.5, 'windflag': 3, 'wdflag': 1, 'alpha1': 1.0, 'pts1': 0.001, 'pts3': 0.02, 'pts2': 0.01, 'epsnov': 0.001, 'hewind': 0.5, 'ck': 1000, 'bwind': 0.0, 'lambdaf': 0.0, 'mxns': 3.0, 'beta': -1.0, 'tflag': 1, 'acc2': 1.5, 'grflag' : 1, 'remnantflag': 4, 'ceflag': 0, 'eddfac': 1.0, 'ifflag': 0, 'bconst': 3000, 'sigma': 265.0, 'gamma': -2.0, 'pisn': 45.0, 'natal_kick_array' : [[-100.0,-100.0,-100.0,-100.0,0.0], [-100.0,-100.0,-100.0,-100.0,0.0]], 'bhsigmafrac' : 1.0, 'polar_kick_angle' : 90, 'qcrit_array' : [0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0], 'cekickflag' : 2, 'cehestarflag' : 0, 'cemergeflag' : 0, 'ecsn' : 2.25, 'ecsn_mlow' : 1.6, 'aic' : 1, 'ussn' : 0, 'sigmadiv' :-20.0, 'qcflag' : 1, 'eddlimflag' : 0, 'fprimc_array' : [2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0], 'bhspinflag' : 0, 'bhspinmag' : 0.0, 'rejuv_fac' : 1.0, 'rejuvflag' : 0, 'htpmb' : 1, 'ST_cr' : 1, 'ST_tide' : 1, 'bdecayfac' : 1, 'rembar_massloss' : 0.5, 'kickflag' : 0, 'zsun' : 0.014, 'bhms_coll_flag' : 0, 'don_lim' : -1, 'acc_lim' : -1, 'rtmsflag' : 0}

    In [6]: bpp, bcm, initC, kick_info = Evolve.evolve(initialbinarytable=single_binary, BSEDict=BSEDict, timestep_conditions =[['kstar_1=14', 'kstar_2<10','dtp=0.1'], ['kstar_2=14', 'kstar_1<10','dtp=0.1']])

    In [16]: print(bcm[['tphys', 'kstar_1', 'kstar_2', 'mass_1', 'mass_2', 'bin_state']])


*******************
restarting a binary
*******************

COSMIC allows you to restart a binary from any point in its evolution from a COSMIC generated bpp array.
Below we provide an example of the same evolutionary track
started from the beginning and three different points in the evolution, once sometime between the beginning and the first object going supernova, once between the first and second supernova, and finally after both supernova

.. code-block:: python

    single_binary = InitialBinaryTable.InitialBinaries(m1=25.543645, m2=20.99784, porb=446.795757, ecc=0.448872, tphysf=13700.0, kstar1=1, kstar2=1, metallicity=0.002)
    BSEDict = {'xi': 1.0, 'bhflag': 1, 'neta': 0.5, 'windflag': 3, 'wdflag': 1, 'alpha1': 1.0, 'pts1': 0.001, 'pts3': 0.02, 'pts2': 0.01, 'epsnov': 0.001, 'hewind': 0.5, 'ck': 1000, 'bwind': 0.0, 'lambdaf': 0.0, 'mxns': 3.0, 'beta': -1.0, 'tflag': 1, 'acc2': 1.5, 'remnantflag': 3, 'ceflag': 0, 'eddfac': 1.0, 'ifflag': 0, 'bconst': 3000, 'sigma': 265.0, 'gamma': -2.0, 'pisn': 45.0, 'natal_kick_array' : [[-100.0,-100.0,-100.0,-100.0,0.0], [-100.0,-100.0,-100.0,-100.0,0.0]], 'bhsigmafrac' : 1.0, 'polar_kick_angle' : 90, 'qcrit_array' : [0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0], 'cekickflag' : 2, 'cehestarflag' : 0, 'cemergeflag' : 0, 'ecsn' : 2.5, 'ecsn_mlow' : 1.4, 'aic' : 1, 'ussn' : 0, 'sigmadiv' :-20.0, 'qcflag' : 1, 'eddlimflag' : 0, 'fprimc_array' : [2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0], 'bhspinflag' : 0, 'bhspinmag' : 0.0, 'rejuv_fac' : 1.0, 'rejuvflag' : 0, 'htpmb' : 1, 'ST_cr' : 1, 'ST_tide' : 0, 'bdecayfac' : 1, 'randomseed' : -1235453, 'grflag' : 1, 'rembar_massloss' : 0.5, 'kickflag' : 0, 'zsun' : 0.014,  'grflag' : 1, 'bhms_coll_flag' : 0, 'don_lim' : -1, 'acc_lim' : -1, 'rtmsflag' : 0}
    for i in [3, 7, 11]:
        bpp, bcm, initC, kick_info = Evolve.evolve(initialbinarytable=single_binary, BSEDict=BSEDict)
        for column in bpp.columns:
            initC = initC.assign(**{column:bpp.iloc[i][column]})
        bpp_mid, bcm_mid, initC_mid, kick_info = Evolve.evolve(initialbinarytable=initC)
        if i == 3:
            print("From beginning")
            print(bpp)
        print("Started in middle at Index {0}".format(i))
        print(bpp_mid)

