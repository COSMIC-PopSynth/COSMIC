.. _timesteps_resolution:

***********************
Dynamic time resolution
***********************

COSMIC has the ability to set time resolution of the bcm array depending on the current state of the evolution. It does this using the ``timestep_conditions`` parameter, which we'll get into below.

How it works
============

Simply put, you tell ``COSMIC`` what time resolution to use during different phases of evolution based on a series of conditions. These conditions can cover everything from the stellar type and masses of the stars, to whether or not they are undergoing mass transfer, to the binary state, etc. Below we'll talk through how to specify those conditions, and then demonstrate their usage in a few worked examples.

Timestep condition syntax
=========================

The timestep resolution for the bcm array can be set using ``timestep_conditions`` parameter. The format of this parameter is a list of up to 15 lists, where each sublist contains zero or more conditions followed by the desired resolution. The conditions are in the form of a string with the format ``'column_name operator condition'``. The resolution is in the form of a string with the format ``'dtp=resolution'``. The columns that can be used are listed in ``cosmic.checkstate.CHECKSTATE_COLUMNS``.

The conditions are evaluated at each time step and if the condition is met, the resolution is set to the desired value. Within a sublist, **all** of the conditions must be met for its resolution to apply. A sublist with no conditions (just a ``dtp``) always applies. If more than one sublist applies at the same time, the **last** one in the list sets the resolution, so put general conditions first and more specific ones after. When no sublist applies, the resolution goes back to the ``dtp`` used for the binary.

The resolution can be any non-negative float value given in units of Myr. If the resolution is set to 0.0, every time step that COSMIC uses internally will be outputted in the bcm array. If the resolution is set to ``None`` then time resolution ends and only the final state is outputted.

Conditions can use the operators ``==``, ``!=``, ``<``, ``<=``, ``>`` and ``>=``. A column can be given both a lower and an upper bound, either as a chained comparison (e.g. ``'1<kstar_1<4'``) or as two separate conditions (e.g. ``'kstar_1>1', 'kstar_1<4'``), as well as one ``!=``. Values can be negative (e.g. ``'sep==-1'``) and spaces around the operators are ignored. For example, the following uses a coarse resolution on the main sequence, a fine resolution on the giant branch, every time step during core helium burning and a medium resolution until the star becomes a remnant

.. code-block:: python

    timestep_conditions = [['kstar_1 < 2', 'dtp=5.0'], ['1 < kstar_1 < 4', 'dtp=0.01'],
                           ['kstar_1 == 4', 'dtp=0.0'], ['4 < kstar_1 < 10', 'dtp=0.5']]

Setting the same bound twice for a column in one sublist (e.g. ``'kstar_1 > 1', 'kstar_1 > 2'``), using an unsupported operator or leaving out the ``dtp`` raises a ``ValueError``.

Syntax examples
---------------

Before we dive in to using these with ``COSMIC``, here are some examples of valid syntax for the ``timestep_conditions`` parameter. Each example describes the condition and then shows the corresponding syntax. Note that the order of the conditions in a sublist does not matter, but the order of the sublists does.

**Basic**

.. code-block:: python

    # increase the timestep resolution during the Hertzsprung gap (``kstar_1==2``)
    timestep_conditions = [['kstar_1 == 2', 'dtp=0.01']]

**Ranged conditions**

.. code-block:: python

    # show all timesteps when a primary star has a specific temperature range
    timestep_conditions = [['10000 < teff_1 < 20000', 'dtp=0.0']]

**Multiple conditions (AND)**

.. code-block:: python

    # show all timesteps when a primary star has a specific temperature range
    # AND it is on its main sequence (``kstar_1<=1``)
    # any of these are equivalent and valid
    timestep_conditions = [['10000 < teff_1 < 20000', 'kstar_1 <= 1', 'dtp=0.0']]
    timestep_conditions = [['10000 < teff_1 < 20000 and kstar_1 <= 1', 'dtp=0.0']]
    timestep_conditions = [['10000 < teff_1 < 20000 & kstar_1 <= 1', 'dtp=0.0']]

**Multiple separate conditions (OR)**

.. code-block:: python

    # use 1 Myr timesteps either when the primary is a WD
    # OR when it's a black hole
    timestep_conditions = [
        ['10 <= kstar_1 <= 12', 'dtp=1.0'],
        ['kstar_2 == 14', 'dtp=1.0']
    ]

**Overlapping conditions** 

.. code-block:: python

    # use 1 Myr timesteps when the primary is a WD,
    # but use 0.1 Myr timesteps when the secondary is a WD
    timestep_conditions = [
        ['10 <= kstar_1 <= 12', 'dtp=1.0'],
        ['10 <= kstar_2 <= 12', 'dtp=0.1']
    ]
    # if *both* stars are WDs, the last condition always takes
    # precedence so the resolution is set to 0.1 Myr

Worked examples with COSMIC
===========================

Below we demonstrate three scenarios, setting dtp only during mass transfer, setting dtp to the same resolution for all of the evolution except for after the system merges or is disrupted, and finally an example of setting dtp only during the HMXB stage of the evolution.

All steps during mass transfer
------------------------------

First, print all time steps during mass transfer

.. ipython:: python

    from cosmic.sample.initialbinarytable import InitialBinaryTable

    from cosmic.evolve import Evolve

.. ipython:: python
    :okwarning:

    single_binary = InitialBinaryTable.InitialBinaries(
        m1=7, m2=5, porb=3000,
        ecc=0.42, tphysf=13700.0, kstar1=1, kstar2=1, metallicity=0.02
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
        m1=80, m2=40, porb=350,
        ecc=0.5, tphysf=13700.0,
        kstar1=1, kstar2=1, metallicity=0.002
    )

    timestep_conditions = [['kstar_1 == 14', 'kstar_2 < 10','dtp=0.1'], ['kstar_2 == 14', 'kstar_1 < 10','dtp=0.1']]

    bpp, bcm, initC, kick_info = Evolve.evolve(
        initialbinarytable=single_binary, 
        BSEDict=BSEDict, SSEDict=SSEDict,
        timestep_conditions=timestep_conditions
    )

    print(bcm[['tphys', 'kstar_1', 'kstar_2', 'mass_1', 'mass_2', 'bin_state']])
