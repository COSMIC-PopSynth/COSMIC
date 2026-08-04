.. _timesteps_modifiers:

****************************
Modifying internal timesteps
****************************

In addition to modifying the resolution of the output data, users can also modify the internal timesteps used by ``COSMIC``.

Changing the default timesteps (not recommended)
------------------------------------------------

By default, ``COSMIC`` timesteps are set by the user-specified settings ``pts1``, ``pts2``, and ``pts3`` (see :ref:`inifile` for more details). These specify timesteps for different evolutionary phases of the binary.

One could edit these directly in the BSEDict as we've seen in other tutorials.

.. ipython:: python

    from cosmic.sample.initialbinarytable import InitialBinaryTable
    from cosmic.evolve import Evolve

    single_binary = InitialBinaryTable.InitialBinaries(
        m1=85.543645, m2=84.99784, porb=446.795757,
        ecc=0.448872, tphysf=13700.0,
        kstar1=1, kstar2=1, metallicity=0.002
    )

    SSEDict = {'stellar_engine': 'sse'}

.. include:: ../../../_generated/default_bsedict.rst

.. ipython:: python

    bpp, bcm, initC, kick_info = Evolve.evolve(
        initialbinarytable=single_binary,
        BSEDict=BSEDict,
        SSEDict=SSEDict
    )

    # take twice as large timesteps for all phases
    BSEDict['pts1'] *= 2
    BSEDict['pts2'] *= 2
    BSEDict['pts3'] *= 2

    bpp_mod, bcm_mod, initC_mod, kick_info_mod = Evolve.evolve(
        initialbinarytable=single_binary,
        BSEDict=BSEDict,
        SSEDict=SSEDict
    )

From inspection of the resulting tables, you may start to see strange behaviour if the timesteps are too large.

Modifying timesteps as a function of mass (recommended)
-------------------------------------------------------

Users can also specify additional modifiers to these timesteps that will be applied during the evolution. The original BSE timesteps can be modified for higher masses, where extrapolation beyond the original stellar grids can lead to unexpected behaviour.

These modifiers are specified as ``dt_mass_modifiers`` in the :meth:`~cosmic.evolve.Evolve.evolve` function when performing evolution. They are written as a list of tuples of the form ``(min, max, modifier)``, where the fractional modifier is applied to the default timesteps for stars with masses between ``min`` and ``max``.

For example, to take half as large timesteps for stars with masses between 50 and 100 solar masses, you can do:

.. ipython:: python

    dt_mass_modifiers = [(50, 100, 0.5)]

    bpp_mod, bcm_mod, initC_mod, kick_info_mod = Evolve.evolve(
        initialbinarytable=single_binary,
        BSEDict=BSEDict,
        SSEDict=SSEDict,
        dt_mass_modifiers=dt_mass_modifiers
    )

The default choice in ``COSMIC`` is to set ``dt_mass_modifiers = [(40, 70, 0.3), (70, np.inf, 0.1)]``, which we find leads to more numerically stable evolution for very massive stars.

Example: BH masses
^^^^^^^^^^^^^^^^^^

Let's go through an example to see how these modifiers can affect the results of our simulations. Specifically, let's explore the relation between initial mass and BH mass for a high-mass grid.

First, we can create a grid of 500 single stars with initial masses between 50 and 150 solar masses.

.. ipython:: python

    import matplotlib.pyplot as plt
    import numpy as np

    N_BINARIES = 500

    # create an initial binary table of single stars
    IBT_bh = InitialBinaryTable.InitialBinaries(
        m1=np.linspace(50, 150, N_BINARIES),
        m2=np.zeros(N_BINARIES),
        porb=np.zeros(N_BINARIES),
        ecc=np.zeros(N_BINARIES),
        kstar1=[1] * N_BINARIES,
        kstar2=[0] * N_BINARIES,
        metallicity=[0.01] * N_BINARIES,
        tphysf=[200] * N_BINARIES,
    )

Now let's copy our BSEDict but turn of pair-instability supernovae (PISN) to simplify the top end of the plot and make sure any differences are just coming from timesteps.

.. ipython:: python

    # switch off PISN to simplify the plot
    no_PISN_BSEDict = BSEDict.copy()
    no_PISN_BSEDict["pisn"] = 0


We can set up 4 different choices of timesteps to compare:

- The default choice in ``COSMIC``: ``dt_mass_modifiers = [(40, 70, 0.3), (70, np.inf, 0.1)]``
- Half the default BSE timesteps: ``dt_mass_modifiers = [(0, np.inf, 0.5)]``
- The default BSE timesteps: ``dt_mass_modifiers = []``
- Double the default BSE timesteps: ``dt_mass_modifiers = [(0, np.inf, 2)]``

.. ipython:: python

    bpp_timesteps = []
    labels = ['Default dt_mass_modifiers', 'Half BSE',
              'BSE defaults', 'Double BSE']
    dt_mass_modifiers_vals = [
        [(40, 70, 0.3), (70, np.inf, 0.1)],
        [(0, np.inf, 0.5)],
        [],
        [(0, np.inf, 2)]
    ]

Now let's use those to evolve our grid of single stars 4 different times and store the resulting BPPs.

.. ipython:: python

    for dtmm in dt_mass_modifiers_vals:
        bpp, _, _, _ = Evolve.evolve(
            initialbinarytable=IBT_bh,
            SSEDict=SSEDict,
            BSEDict=BSEDict,
            dt_mass_modifiers=dtmm,
        )
        bpp_timesteps.append(bpp)

With our evolution done, let's find the final BH mass and the corresponding initial mass in each grid for plotting.

.. ipython:: python

    m_inits = []
    m_bhs = []

    for bpp in bpp_timesteps:
        bh = bpp.loc[bpp['kstar_1'] == 14].groupby(level=0).last()
        init = IBT_bh.loc[bh.index, 'mass_1']
        m_init, m_bh = init.values, bh['mass_1'].values
        m_inits.append(m_init)
        m_bhs.append(m_bh)

And finally, let's plot the results!

.. ipython:: python

    plt.rc('font', family='serif')
    plt.rcParams['text.usetex'] = False
    fs = 24
    params = {'figure.figsize': (12, 8),
            'legend.fontsize': 0.7*fs,
            'legend.title_fontsize': 0.8*fs,
            'axes.labelsize': fs,
            'xtick.labelsize': 0.9 * fs,
            'ytick.labelsize': 0.9 * fs,
            'axes.linewidth': 1.1,
            'xtick.major.size': 7,
            'xtick.minor.size': 4,
            'ytick.major.size': 7,
            'ytick.minor.size': 4}
    plt.rcParams.update(params)

    fig, axes = plt.subplots(2, 1, figsize=(14, 8), gridspec_kw={'height_ratios': [3, 1]})

    for m_init, m_bh, label, zorder in zip(m_inits, m_bhs, labels, [100, 50, 25, 1]):
        axes[0].plot(m_init, m_bh, zorder=zorder, label=label)
    axes[0].legend();

    for m_init, m_bh in zip(m_inits, m_bhs):
        axes[1].plot(m_init, m_bh - m_bhs[0])
    axes[1].set_xlabel('Initial primary mass $M_1^\\mathrm{init}$ [$M_\\odot$]');
    
    for ax in axes:
        ax.grid(True, lw=0.4, alpha=0.4);
    axes[0].set_ylabel('BH mass [$M_\\odot$]');
    axes[1].set_ylabel('Difference [$M_\\odot$]');

    plt.tight_layout();
    
    @savefig dt_mass_modifiers.png
    plt.show()

So we can see here that there are clear numerical artifacts in the relation between initial mass and BH mass when using the default BSE timesteps (and especially when they are larger!). This motivated our choice of the default ``dt_mass_modifiers`` in ``COSMIC`` to ensure more stable evolution for very massive stars.