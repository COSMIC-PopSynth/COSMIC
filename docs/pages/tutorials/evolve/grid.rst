*********************************
Evolving a fixed grid of binaries
*********************************

Sometimes it is helpful to run a fixed grid of initial binaries to explore how
changing a single parameter affects the evolved binary. Let's start by importing the necessary modules and
setting up the BSEDict settings as we've done in the previous examples.

.. ipython:: python

    from cosmic.sample.initialbinarytable import InitialBinaryTable

    from cosmic.evolve import Evolve

.. include:: ../../../_generated/default_bsedict.rst


Here we evolve the same system that produces a GW150914-like binary, but run over several initial orbital
periods spaced evenly in log space.

.. ipython:: python
    :okwarning:

    n_grid = 10

    binary_grid = InitialBinaryTable.InitialBinaries(
        m1=np.ones(n_grid)*100.0,
        m2=np.ones(n_grid)*85.0,
        porb=np.logspace(3,5,n_grid),
        ecc=np.ones(n_grid)*0.65,
        tphysf=np.ones(n_grid)*13700.0,
        kstar1=np.ones(n_grid),
        kstar2=np.ones(n_grid),
        metallicity=np.ones(n_grid)*0.005
    )

    SSEDict = {'stellar_engine': 'sse'}

    print(binary_grid)

    bpp, bcm, initC, kick_info = Evolve.evolve(initialbinarytable=binary_grid, BSEDict=BSEDict, SSEDict=SSEDict)

    print(bpp)

    print(bcm)
