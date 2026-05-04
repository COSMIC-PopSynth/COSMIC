************************
Evolving a single binary
************************

Initial conditions
==================

Below is the process to initialize and evolve a binary that
could have formed a GW150914-like binary. First, import the modules in COSMIC
that initialize and evolve the binary.

.. ipython:: python

    from cosmic.sample.initialbinarytable import InitialBinaryTable
    from cosmic.evolve import Evolve


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

.. ipython:: python

    single_binary = InitialBinaryTable.InitialBinaries(m1=85.543645, m2=84.99784, porb=446.795757,
                                                    ecc=0.448872, tphysf=13700.0,
                                                    kstar1=1, kstar2=1, metallicity=0.002)

    print(single_binary)


There are two available methods for evaluating the evolution of each star in the binary:

* SSE - the single star evolution fitting formulae from `Hurley+2000 <https://ui.adsabs.harvard.edu/abs/2000MNRAS.315..543H/abstract>`_

* METISSE - a new stellar evolution package that uses interpolation of detailed stellar models supplied by the user. 

You can learn more these choices by consulting the `METISSE documentation <https://metisse.readthedocs.io/en/latest/>`_.

If you have a grid of detailed stellar models and want to use METISSE, you can specify the
SSEDict as follows: 

.. ipython::

    In [3]: SSEDict = {'stellar_engine': 'metisse',
       ...:             'path_to_tracks': 'path/to/your/hydrogen/rich/stellar/models/',
       ...:             'path_to_he_tracks': 'path/to/your/helium/stellar/models/'}


For now, we will use the default SSE method for evolving single stars in binaries. 
To do this, we specify the SSEDict which tells COSMIC to use the Hurley fitting formulae.

.. ipython::

    In [3]: SSEDict = {'stellar_engine': 'sse'} 

(Binary) stellar physics assumptions
====================================

The flags for the various binary evolution prescriptions used in BSE also need to be set.
Each flag is saved in the BSEDict dictionary. Note that the BSEDict
only needs to be specified the first time a binary is evolved with COSMIC or
if you need to change the binary evolution prescriptions.

If you are unfamiliar with these prescriptions, it is 
advised to run the defaults from the COSMIC install which are consistent
with `Breivik+2020 <https://ui.adsabs.harvard.edu/abs/2019arXiv191100903B/abstract>`_,
though we don't promise that these are the most up-to-date prescriptions.

.. include:: ../../_generated/default_bsedict.rst


Running a binary
================

Once the binary is initialized and the BSE model is set, the system is evolved with the
the Evolve class, which calls the evolv2.f subroutine in the BSE source code.

.. ipython:: python
    :okwarning:

    bpp, bcm, initC, kick_info = Evolve.evolve(initialbinarytable=single_binary, BSEDict=BSEDict, SSEDict=SSEDict)


Output
======

For every evolved binary system, BSE generates two arrays, which are stored as pandas DataFrames in COSMIC:

* bpp - contains binary parameters at important stages in the binary's evolution, including stellar evolutionary phase changes or mass transfer episodes.

* bcm - contains several binary parameters at user specified time steps during the binary's evolution. The default setting in COSMIC is to output the final stage of the binary at the evolution time specified by the user.

You can see the different parameters included in each DataFrame using the columns attribute of the DataFrame:

.. ipython:: python
    
    print(bpp.columns)

    print(bcm.columns)

The units are broadly consistent with BSE and are described in :ref:`output_info`.

The evol_type column in bpp indicates the evolutionary change that occurred for each line.
The meaning of each number is described here, :ref:`evolve-type-table`.

Each of the parameters in bpp or bcm can be accessed in the usual way for DataFrames:

.. ipython:: python

    print(bpp.mass_1)

    print(bpp[['mass_1', 'mass_2', 'kstar_1', 'kstar_2', 'sep', 'evol_type']])


You can use the ``utils.convert_kstar_evol_type`` function to convert the
``kstar_1``, ``kstar_2``, and ``evol_type`` columns from integers to strings
that describe each int:

.. ipython:: python

    from cosmic.utils import convert_kstar_evol_type

    convert_kstar_evol_type(bpp)


Note that ``utils.convert_kstar_evol_type`` is only applicable to the bpp
array.

Modifying the columns in each table
-----------------------------------

The columns in each table can be modified by passing in a list of desired columns to the
``bpp_columns`` or ``bcm_columns`` keyword arguments in the ``Evolve.evolve`` method.
This is useful if you only want a subset of the available columns to reduce memory usage, or if you
want columns from ``bcm`` in the ``bpp`` table or vice versa.

For example, to only get the time, masses, stellar types, separation, and evolution type, you can do:

.. ipython:: python
    
    bpp, bcm, initC, kick_info = Evolve.evolve(
        initialbinarytable=single_binary,
        BSEDict=BSEDict,
        bpp_columns=['tphys', 'mass_1', 'mass_2', 'kstar_1', 'kstar_2', 'sep', 'evol_type']
    )

    print(bpp)

.. note::

    Whichever columns are specified in ``bpp_columns`` or ``bcm_columns``, there will always be a ``bin_num``
    column included in the output tables to identify each binary uniquely.

Plotting the evolution
======================

You can also use the built-in plotting function to see how the system evolves:

.. ipython:: python
    :okwarning:

    from cosmic.plotting import evolve_and_plot
    fig = evolve_and_plot(single_binary, t_min=None, t_max=None, BSEDict=BSEDict, SSEDict=SSEDict, sys_obs={})

.. plot::

    from cosmic.sample.initialbinarytable import InitialBinaryTable
    from cosmic.plotting import evolve_and_plot
    single_binary = InitialBinaryTable.InitialBinaries(m1=85.543645, m2=84.99784, porb=446.795757, ecc=0.448872, tphysf=13700.0, kstar1=1, kstar2=1, metallicity=0.002)
    SSEDict = {'stellar_engine': 'sse'}
    fig = evolve_and_plot(single_binary, t_min=None, t_max=None, BSEDict=default_BSEDict, SSEDict=SSEDict, sys_obs={})


In this case, all the action happens in the first few Myr, so let's specify a t_max:

.. ipython:: python
    :okwarning:
    
    fig = evolve_and_plot(initC, t_min=None, t_max=6.0, BSEDict=BSEDict, SSEDict=SSEDict, sys_obs={})

.. plot::

    from cosmic.sample.initialbinarytable import InitialBinaryTable
    from cosmic.plotting import evolve_and_plot
    single_binary = InitialBinaryTable.InitialBinaries(m1=85.543645, m2=84.99784, porb=446.795757, ecc=0.448872, tphysf=13700.0, kstar1=1, kstar2=1, metallicity=0.002)
    SSEDict = {'stellar_engine': 'sse'}
    fig = evolve_and_plot(single_binary, t_min=None, t_max=6.0, BSEDict=default_BSEDict, SSEDict=SSEDict, sys_obs={})
