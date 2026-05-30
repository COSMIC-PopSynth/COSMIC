**************************
Evolving multiple binaries
**************************

Following on from `evolving a single binary <single.html>`_, we can also evolve multiple binaries.
Let's start by importing the necessary modules:

.. ipython:: python

    from cosmic.sample.initialbinarytable import InitialBinaryTable

    from cosmic.evolve import Evolve


And use the same SSE and BSE dictionaries as before:

.. ipython:: python

    SSEDict = {'stellar_engine': 'sse'}

.. include:: ../../../_generated/default_bsedict.rst

Below is an example for systems that could form GW150914 and GW170817 - like binaries.

.. ipython:: python
    :okwarning:

    binary_set = InitialBinaryTable.InitialBinaries(
        m1=[85.543645, 11.171469], m2=[84.99784, 6.67305],
        porb=[446.795757, 170.758343], ecc=[0.448872, 0.370],
        tphysf=[13700.0, 13700.0],
        kstar1=[1, 1], kstar2=[1, 1],
        metallicity=[0.002, 0.02]
    )

    print(binary_set)

    import numpy as np

    np.random.seed(5)

    bpp, bcm, initC, kick_info = Evolve.evolve(initialbinarytable=binary_set, BSEDict=BSEDict, SSEDict=SSEDict)

As before, bpp, bcm, and initC are returned as pandas DataFrames which assign an
index to each binary system we evolve. We can access each binary as follows:

.. ipython:: python

    print(bpp.loc[0])

    print(bcm.loc[0])

    print(initC.loc[0])

    print(bpp.loc[1])

The plotting function can also take in multiple binaries. Let's plot both the GW150914-like
progenitor evolution and the GW170817-like progenitor evolutions. For the GW170817-like
progenitor, we expect most of the evolution to take place in the first ~60 Myr.

.. ipython:: python
    :okwarning:
    :okexcept:

    from cosmic.plotting import evolve_and_plot
    fig = evolve_and_plot(binary_set, t_min=None, t_max=[6.0, 60.0], BSEDict=BSEDict, SSEDict=SSEDict, sys_obs={})


.. plot::

    from cosmic.sample.initialbinarytable import InitialBinaryTable
    from cosmic.plotting import evolve_and_plot
    import numpy as np
    np.random.seed(5)
    binary_set = InitialBinaryTable.InitialBinaries(m1=[85.543645, 11.171469], m2=[84.99784, 9.67305], porb=[446.795757, 370.758343], ecc=[0.448872, 0.370], tphysf=[13700.0, 13700.0], kstar1=[1, 1], kstar2=[1, 1], metallicity=[0.002, 0.02])
    fig = evolve_and_plot(binary_set, t_min=None, t_max=[6.0, 60.0], BSEDict=default_BSEDict, SSEDict=default_SSEDict, sys_obs={})

