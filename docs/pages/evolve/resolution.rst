***********************
Dynamic time resolution
***********************

COSMIC has the ability to set time resolution of the bcm array depending on the current state of the evolution.

.. ipython::

    In [1]: from cosmic.sample.initialbinarytable import InitialBinaryTable

    In [2]: from cosmic.evolve import Evolve

General syntax
==============

The timestep resolution for the bcm array can be set using ``timestep_conditions`` parameter. The format of
this parameter is a list of lists, where each sublist contains zero or more conditions followed by the desired
resolution. The conditions are in the form of a string with the format ``'column_name=condition'``. The resolution
is in the form of a string with the format ``'dtp=resolution'``.

The conditions are evaluated at each time step and if the condition is met, the resolution is set to the desired value. If multiple conditions are specified, a timestep is outputted if **any** of the conditions are met.

The resolution can be any positive float value given in units of Myr. If the resolution is set to 0.0, every time step will be outputted in the bcm array.

Examples
========

Below we demonstrate three scenarios, setting dtp only during mass transfer, setting dtp to the same resolution for all of the evolution except for after the system merges or is disrupted, and finally an example of setting dtp only during the HMXB stage of the evolution.

All steps during mass transfer
------------------------------

First, print all time steps during mass transfer

.. ipython::
    :okwarning:

    In [3]: single_binary = InitialBinaryTable.InitialBinaries(m1=7.806106, m2=5.381412, porb=2858.942021,
       ...:                                                    ecc=0.601408, tphysf=13700.0, kstar1=1, kstar2=1, metallicity=0.02)

    In [4]: BSEDict = {'xi': 1.0, 'bhflag': 1, 'neta': 0.5, 'windflag': 3, 'wdflag': 1, 'alpha1': 1.0, 'pts1': 0.001, 'pts3': 0.02, 'pts2': 0.01, 'epsnov': 0.001, 'hewind': 0.5, 'ck': 1000, 'bwind': 0.0, 'lambdaf': 0.0, 'mxns': 3.0, 'beta': -1.0, 'tflag': 1, 'acc2': 1.5, 'grflag' : 1, 'remnantflag': 4, 'ceflag': 0, 'eddfac': 1.0, 'ifflag': 0, 'bconst': 3000, 'sigma': 265.0, 'gamma': -2.0, 'pisn': 45.0, 'natal_kick_array' : [[-100.0,-100.0,-100.0,-100.0,0.0], [-100.0,-100.0,-100.0,-100.0,0.0]], 'bhsigmafrac' : 1.0, 'polar_kick_angle' : 90, 'qcrit_array' : [0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0], 'cekickflag' : 2, 'cehestarflag' : 0, 'cemergeflag' : 0, 'ecsn' : 2.25, 'ecsn_mlow' : 1.6, 'aic' : 1, 'ussn' : 0, 'sigmadiv' :-20.0, 'qcflag' : 1, 'eddlimflag' : 0, 'fprimc_array' : [2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0], 'bhspinflag' : 0, 'bhspinmag' : 0.0, 'rejuv_fac' : 1.0, 'rejuvflag' : 0, 'htpmb' : 1, 'ST_cr' : 1, 'ST_tide' : 1, 'bdecayfac' : 1, 'rembar_massloss' : 0.5, 'kickflag' : 1, 'zsun' : 0.014, 'bhms_coll_flag' : 0, 'don_lim' : -1, 'acc_lim' : -1, 'rtmsflag' : 0, 'wd_mass_lim': 1}

    In [5]: # Define the condition for the time step
       ...: timestep_conditions = [['RRLO_1>=1', 'dtp=0.0'], ['RRLO_2>=1', 'dtp=0.0']]

    In [6]: bpp, bcm, initC, kick_info = Evolve.evolve(initialbinarytable=single_binary, BSEDict=BSEDict,
       ...:                                            timestep_conditions=timestep_conditions)

    In [7]: print(bcm[['tphys', 'kstar_1', 'kstar_2', 'mass_1', 'mass_2', 'RRLO_1', 'RRLO_2']])


Condition on binary state
-------------------------
Second, pick a certain resolution for the bcm array until the system merges or is disrupted and then only print the final state

.. ipython::
    :okwarning:

    In [8]: timestep_conditions = [['binstate=0', 'dtp=1.0']]

    In [9]: bpp, bcm, initC, kick_info = Evolve.evolve(initialbinarytable=single_binary, BSEDict=BSEDict,
       ...:                                            timestep_conditions=timestep_conditions)

    In [10]: print(bcm[['tphys', 'kstar_1', 'kstar_2', 'mass_1', 'mass_2', 'bin_state']])


Condition on evolutionary state
-------------------------------
Finally, we show how to print a fine resolution only during the HMXB stage of the evolution.

.. ipython::
    :okwarning:

    In [11]: single_binary = InitialBinaryTable.InitialBinaries(m1=85.543645, m2=84.99784, porb=446.795757,
       ....:                                                    ecc=0.448872, tphysf=13700.0,
       ....:                                                    kstar1=1, kstar2=1, metallicity=0.002)

    In [12]: BSEDict = {'xi': 1.0, 'bhflag': 1, 'neta': 0.5, 'windflag': 3, 'wdflag': 1, 'alpha1': 1.0, 'pts1': 0.001, 'pts3': 0.02, 'pts2': 0.01, 'epsnov': 0.001, 'hewind': 0.5, 'ck': 1000, 'bwind': 0.0, 'lambdaf': 0.0, 'mxns': 3.0, 'beta': -1.0, 'tflag': 1, 'acc2': 1.5, 'grflag' : 1, 'remnantflag': 4, 'ceflag': 0, 'eddfac': 1.0, 'ifflag': 0, 'bconst': 3000, 'sigma': 265.0, 'gamma': -2.0, 'pisn': 45.0, 'natal_kick_array' : [[-100.0,-100.0,-100.0,-100.0,0.0], [-100.0,-100.0,-100.0,-100.0,0.0]], 'bhsigmafrac' : 1.0, 'polar_kick_angle' : 90, 'qcrit_array' : [0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0], 'cekickflag' : 2, 'cehestarflag' : 0, 'cemergeflag' : 0, 'ecsn' : 2.25, 'ecsn_mlow' : 1.6, 'aic' : 1, 'ussn' : 0, 'sigmadiv' :-20.0, 'qcflag' : 1, 'eddlimflag' : 0, 'fprimc_array' : [2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0], 'bhspinflag' : 0, 'bhspinmag' : 0.0, 'rejuv_fac' : 1.0, 'rejuvflag' : 0, 'htpmb' : 1, 'ST_cr' : 1, 'ST_tide' : 1, 'bdecayfac' : 1, 'rembar_massloss' : 0.5, 'kickflag' : 1, 'zsun' : 0.014, 'bhms_coll_flag' : 0, 'don_lim' : -1, 'acc_lim' : -1, 'rtmsflag' : 0, 'wd_mass_lim': 1}

    In [13]: timestep_conditions = [['kstar_1=14', 'kstar_2<10','dtp=0.1'], ['kstar_2=14', 'kstar_1<10','dtp=0.1']]

    In [14]: bpp, bcm, initC, kick_info = Evolve.evolve(initialbinarytable=single_binary, BSEDict=BSEDict,
       ....:                                            timestep_conditions=timestep_conditions)

    In [15]: print(bcm[['tphys', 'kstar_1', 'kstar_2', 'mass_1', 'mass_2', 'bin_state']])
