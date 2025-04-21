************************
Evolving a single binary
************************

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

    In [3]: single_binary = InitialBinaryTable.InitialBinaries(m1=85.543645, m2=84.99784, porb=446.795757,
       ...:                                                    ecc=0.448872, tphysf=13700.0,
       ...:                                                    kstar1=1, kstar2=1, metallicity=0.002)

    In [4]: print(single_binary)


The flags for the various binary evolution prescriptions used in BSE also need to be set.
Each flag is saved in the BSEDict dictionary. Note that the BSEDict
only needs to be specified the first time a binary is evolved with COSMIC or
if you need to change the binary evolution prescriptions.

If you are unfamiliar with these prescriptions, it is highly
advised to either run the defaults from the COSMIC install which are consistent
with `Breivik+2020 <https://ui.adsabs.harvard.edu/abs/2019arXiv191100903B/abstract>`_

.. ipython::

    In [5]: BSEDict = {'xi': 1.0, 'bhflag': 1, 'neta': 0.5, 'windflag': 3, 'wdflag': 1, 'alpha1': 1.0, 'pts1': 0.001, 'pts3': 0.02, 'pts2': 0.01, 'epsnov': 0.001, 'hewind': 0.5, 'ck': 1000, 'bwind': 0.0, 'lambdaf': 0.0, 'mxns': 3.0, 'beta': -1.0, 'tflag': 1, 'acc2': 1.5, 'grflag' : 1, 'remnantflag': 4, 'ceflag': 0, 'eddfac': 1.0, 'ifflag': 0, 'bconst': 3000, 'sigma': 265.0, 'gamma': -2.0, 'pisn': 45.0, 'natal_kick_array' : [[-100.0,-100.0,-100.0,-100.0,0.0], [-100.0,-100.0,-100.0,-100.0,0.0]], 'bhsigmafrac' : 1.0, 'polar_kick_angle' : 90, 'qcrit_array' : [0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0], 'cekickflag' : 2, 'cehestarflag' : 0, 'cemergeflag' : 0, 'ecsn' : 2.25, 'ecsn_mlow' : 1.6, 'aic' : 1, 'ussn' : 0, 'sigmadiv' :-20.0, 'qcflag' : 1, 'eddlimflag' : 0, 'fprimc_array' : [2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0], 'bhspinflag' : 0, 'bhspinmag' : 0.0, 'rejuv_fac' : 1.0, 'rejuvflag' : 0, 'htpmb' : 1, 'ST_cr' : 1, 'ST_tide' : 1, 'bdecayfac' : 1, 'rembar_massloss' : 0.5, 'kickflag' : 1, 'zsun' : 0.014, 'bhms_coll_flag' : 0, 'don_lim' : -1, 'acc_lim' : -1, 'rtmsflag' : 0, 'wd_mass_lim': 1}


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

You can also use the built-in plotting function to see how the system evolves:

.. ipython::
    :okwarning:

    In [13]: from cosmic.plotting import evolve_and_plot

    In [14]: single_binary = InitialBinaryTable.InitialBinaries(m1=85.543645, m2=84.99784, porb=446.795757,
       ....:                                                    ecc=0.448872, tphysf=13700.0,
       ....:                                                    kstar1=1, kstar2=1, metallicity=0.002)

    In [15]: BSEDict = {'xi': 1.0, 'bhflag': 1, 'neta': 0.5, 'windflag': 3, 'wdflag': 1, 'alpha1': 1.0, 'pts1': 0.001, 'pts3': 0.02, 'pts2': 0.01, 'epsnov': 0.001, 'hewind': 0.5, 'ck': 1000, 'bwind': 0.0, 'lambdaf': 0.0, 'mxns': 3.0, 'beta': -1.0, 'tflag': 1, 'acc2': 1.5, 'grflag' : 1, 'remnantflag': 4, 'ceflag': 0, 'eddfac': 1.0, 'ifflag': 0, 'bconst': 3000, 'sigma': 265.0, 'gamma': -2.0, 'pisn': 45.0, 'natal_kick_array' : [[-100.0,-100.0,-100.0,-100.0,0.0], [-100.0,-100.0,-100.0,-100.0,0.0]], 'bhsigmafrac' : 1.0, 'polar_kick_angle' : 90, 'qcrit_array' : [0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0], 'cekickflag' : 2, 'cehestarflag' : 0, 'cemergeflag' : 0, 'ecsn' : 2.25, 'ecsn_mlow' : 1.6, 'aic' : 1, 'ussn' : 0, 'sigmadiv' :-20.0, 'qcflag' : 1, 'eddlimflag' : 0, 'fprimc_array' : [2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0], 'bhspinflag' : 0, 'bhspinmag' : 0.0, 'rejuv_fac' : 1.0, 'rejuvflag' : 0, 'htpmb' : 1, 'ST_cr' : 1, 'ST_tide' : 1, 'bdecayfac' : 1, 'rembar_massloss' : 0.5, 'kickflag' : 1, 'zsun' : 0.014, 'bhms_coll_flag' : 0, 'don_lim' : -1, 'acc_lim' : -1, 'rtmsflag' : 0, 'wd_mass_lim': 1}

    In [16]: fig = evolve_and_plot(single_binary, t_min=None, t_max=None, BSEDict=BSEDict, sys_obs={})

.. plot::

    from cosmic.sample.initialbinarytable import InitialBinaryTable
    from cosmic.plotting import evolve_and_plot
    single_binary = InitialBinaryTable.InitialBinaries(m1=85.543645, m2=84.99784, porb=446.795757, ecc=0.448872, tphysf=13700.0, kstar1=1, kstar2=1, metallicity=0.002)
    BSEDict = {'xi': 1.0, 'bhflag': 1, 'neta': 0.5, 'windflag': 3, 'wdflag': 1, 'alpha1': 1.0, 'pts1': 0.001, 'pts3': 0.02, 'pts2': 0.01, 'epsnov': 0.001, 'hewind': 0.5, 'ck': 1000, 'bwind': 0.0, 'lambdaf': 0.0, 'mxns': 3.0, 'beta': -1.0, 'tflag': 1, 'acc2': 1.5, 'grflag' : 1, 'remnantflag': 4, 'ceflag': 0, 'eddfac': 1.0, 'ifflag': 0, 'bconst': 3000, 'sigma': 265.0, 'gamma': -2.0, 'pisn': 45.0, 'natal_kick_array' : [[-100.0,-100.0,-100.0,-100.0,0.0], [-100.0,-100.0,-100.0,-100.0,0.0]], 'bhsigmafrac' : 1.0, 'polar_kick_angle' : 90, 'qcrit_array' : [0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0], 'cekickflag' : 2, 'cehestarflag' : 0, 'cemergeflag' : 0, 'ecsn' : 2.25, 'ecsn_mlow' : 1.6, 'aic' : 1, 'ussn' : 0, 'sigmadiv' :-20.0, 'qcflag' : 1, 'eddlimflag' : 0, 'fprimc_array' : [2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0], 'bhspinflag' : 0, 'bhspinmag' : 0.0, 'rejuv_fac' : 1.0, 'rejuvflag' : 0, 'htpmb' : 1, 'ST_cr' : 1, 'ST_tide' : 1, 'bdecayfac' : 1, 'rembar_massloss' : 0.5, 'kickflag' : 1, 'zsun' : 0.014, 'bhms_coll_flag' : 0, 'don_lim' : -1, 'acc_lim' : -1, 'rtmsflag' : 0, 'wd_mass_lim': 1}
    fig = evolve_and_plot(single_binary, t_min=None, t_max=None, BSEDict=BSEDict, sys_obs={})


In this case, all the action happens in the first few Myr, so let's specify a t_max:

.. ipython::
    :okwarning:

    In [17]: fig = evolve_and_plot(initC, t_min=None, t_max=6.0, BSEDict={}, sys_obs={})

.. plot::

    from cosmic.sample.initialbinarytable import InitialBinaryTable
    from cosmic.plotting import evolve_and_plot
    single_binary = InitialBinaryTable.InitialBinaries(m1=85.543645, m2=84.99784, porb=446.795757, ecc=0.448872, tphysf=13700.0, kstar1=1, kstar2=1, metallicity=0.002)
    BSEDict = {'xi': 1.0, 'bhflag': 1, 'neta': 0.5, 'windflag': 3, 'wdflag': 1, 'alpha1': 1.0, 'pts1': 0.001, 'pts3': 0.02, 'pts2': 0.01, 'epsnov': 0.001, 'hewind': 0.5, 'ck': 1000, 'bwind': 0.0, 'lambdaf': 0.0, 'mxns': 3.0, 'beta': -1.0, 'tflag': 1, 'acc2': 1.5, 'grflag' : 1, 'remnantflag': 4, 'ceflag': 0, 'eddfac': 1.0, 'ifflag': 0, 'bconst': 3000, 'sigma': 265.0, 'gamma': -2.0, 'pisn': 45.0, 'natal_kick_array' : [[-100.0,-100.0,-100.0,-100.0,0.0], [-100.0,-100.0,-100.0,-100.0,0.0]], 'bhsigmafrac' : 1.0, 'polar_kick_angle' : 90, 'qcrit_array' : [0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0], 'cekickflag' : 2, 'cehestarflag' : 0, 'cemergeflag' : 0, 'ecsn' : 2.25, 'ecsn_mlow' : 1.6, 'aic' : 1, 'ussn' : 0, 'sigmadiv' :-20.0, 'qcflag' : 1, 'eddlimflag' : 0, 'fprimc_array' : [2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0], 'bhspinflag' : 0, 'bhspinmag' : 0.0, 'rejuv_fac' : 1.0, 'rejuvflag' : 0, 'htpmb' : 1, 'ST_cr' : 1, 'ST_tide' : 1, 'bdecayfac' : 1, 'rembar_massloss' : 0.5, 'kickflag' : 1, 'zsun' : 0.014, 'bhms_coll_flag' : 0, 'don_lim' : -1, 'acc_lim' : -1, 'rtmsflag' : 0, 'wd_mass_lim': 1}
    fig = evolve_and_plot(single_binary, t_min=None, t_max=6.0, BSEDict=BSEDict, sys_obs={})