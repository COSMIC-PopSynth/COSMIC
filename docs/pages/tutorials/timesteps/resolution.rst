.. _timesteps_resolution:

***********************
Dynamic time resolution
***********************

COSMIC has the ability to set time resolution of the bcm array depending on the current state of the evolution.

.. ipython:: python

    from cosmic.sample.initialbinarytable import InitialBinaryTable

    from cosmic.evolve import Evolve

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

.. ipython:: python
    :okwarning:

    single_binary = InitialBinaryTable.InitialBinaries(
        m1=7.806106, m2=5.381412, porb=2858.942021,
        ecc=0.601408, tphysf=13700.0, kstar1=1, kstar2=1, metallicity=0.02
    )

    SSEDict = {'stellar_engine': 'sse'}

.. include:: ../../../_generated/default_bsedict.rst

.. ipython:: python
    :okwarning:

    # Define the condition for the time step
    timestep_conditions = [['RRLO_1>=1', 'dtp=0.0'], ['RRLO_2>=1', 'dtp=0.0']]

    bpp, bcm, initC, kick_info = Evolve.evolve(
        initialbinarytable=single_binary, BSEDict=BSEDict, SSEDict=SSEDict,
        timestep_conditions=timestep_conditions
    )
    print(bcm[['tphys', 'kstar_1', 'kstar_2', 'mass_1', 'mass_2', 'RRLO_1', 'RRLO_2']])


Condition on binary state
-------------------------
Second, pick a certain resolution for the bcm array until the system merges or is disrupted and then only print the final state

.. ipython:: python
    :okwarning:

    timestep_conditions = [['binstate=0', 'dtp=1.0']]

    bpp, bcm, initC, kick_info = Evolve.evolve(
        initialbinarytable=single_binary, 
        BSEDict=BSEDict, SSEDict=SSEDict, 
        timestep_conditions=timestep_conditions
    )

    print(bcm[['tphys', 'kstar_1', 'kstar_2', 'mass_1', 'mass_2', 'bin_state']])


Condition on evolutionary state
-------------------------------
Finally, we show how to print a fine resolution only during the HMXB stage of the evolution.

.. ipython:: python
    :okwarning:

    single_binary = InitialBinaryTable.InitialBinaries(
        m1=85.543645, m2=84.99784, porb=446.795757,
        ecc=0.448872, tphysf=13700.0,
        kstar1=1, kstar2=1, metallicity=0.002
    )

    timestep_conditions = [['kstar_1=14', 'kstar_2<10','dtp=0.1'], ['kstar_2=14', 'kstar_1<10','dtp=0.1']]

    bpp, bcm, initC, kick_info = Evolve.evolve(
        initialbinarytable=single_binary, 
        BSEDict=BSEDict, SSEDict=SSEDict,
        timestep_conditions=timestep_conditions
    )

    print(bcm[['tphys', 'kstar_1', 'kstar_2', 'mass_1', 'mass_2', 'bin_state']])
