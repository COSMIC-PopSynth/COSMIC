.. _rerun_restart:

*******************
Restarting a binary
*******************

COSMIC allows you to restart a binary from any point in its evolution from a COSMIC generated bpp array.

.. ipython:: python

    from cosmic.sample.initialbinarytable import InitialBinaryTable

    from cosmic.evolve import Evolve

Below we provide an example of the same evolutionary track
started from the beginning and three different points in the evolution:

- sometime between the beginning and the first object going supernova
- between the first and second supernova
- after both supernova

.. ipython:: python
    :okwarning:

    single_binary = InitialBinaryTable.InitialBinaries(
        m1=25, m2=20, porb=6000, ecc=0.0,
        tphysf=13700.0, kstar1=1, kstar2=1, metallicity=0.002
    )
    
    SSEDict = {'stellar_engine': 'sse'}
    
.. include:: ../../../_generated/default_bsedict.rst

.. ipython:: python
    :okwarning:
    
    # evolve the binary
    bpp, bcm, initC, kick_info = Evolve.evolve(initialbinarytable=single_binary, BSEDict=BSEDict, SSEDict=SSEDict)
    print("From beginning")
    print(bpp)
    
    # assign row indices and find SNe
    bpp["row_iloc"] = list(range(len(bpp)))
    first_SN_iloc = bpp[bpp["evol_type"] == 15].iloc[0]["row_iloc"]
    second_SN_iloc = bpp[bpp["evol_type"] == 16].iloc[0]["row_iloc"]
    print("First SN Index: ", first_SN_iloc)
    print("Second SN Index: ", second_SN_iloc)

    # choose some inds to restart from
    before_SN_1 = int(first_SN_iloc // 2)
    between_SNe = int((first_SN_iloc + second_SN_iloc) // 2)
    after_SN_2 = int((second_SN_iloc + bpp["row_iloc"].max()) // 2)
    print("Restart Indices: ", before_SN_1, between_SNe, after_SN_2)

.. ipython:: python
    :okwarning:

    # restart from different points
    for i in [before_SN_1, between_SNe, after_SN_2]:
        new_initC = initC.copy()
        for column in bpp.columns:
            new_initC = new_initC.assign(**{column:bpp.iloc[i][column]})
        bpp_mid, bcm_mid, initC_mid, kick_info = Evolve.evolve(initialbinarytable=new_initC, BSEDict={})
        print("Started in middle at Index {0}".format(i))
        print(bpp_mid)

Natal kick example
==================

One example of where restarting a binary can be extremely helpful is studying
how natal kicks affect a binary independently of its previous evolution. This is
particularly relevant for
`Gaia BH1 <https://ui.adsabs.harvard.edu/abs/2023MNRAS.518.1057E/abstract>`_
and `Gaia BH2 <https://ui.adsabs.harvard.edu/abs/2023MNRAS.521.4323E/abstract>`_ which are difficult to produce through the standard common envelope
channels. We can still study the effect of natal kicks on these binaries if we
restart the evolution after the mass transfer would occur. We can do this by using a binary which gets us to the right masses given the metallicity, then overwrite some of the initial conditions to resample the natal kicks and pre-explosion separation. 

.. ipython:: python
    :okwarning:

    from cosmic import utils
    import pandas as pd

    single_binary = InitialBinaryTable.InitialBinaries(
        m1=65.0, m2=0.93, porb=4500, ecc=0.448872, 
        tphysf=13700.0, kstar1=1, kstar2=1, metallicity=0.014*0.6
    )
    bpp, bcm, initC, kick_info = Evolve.evolve(
        initialbinarytable=single_binary, BSEDict=BSEDict, SSEDict=SSEDict
    )

    for column in bpp.columns:
        initC = initC.assign(**{column:bpp.iloc[6][column]})

    initC = pd.concat([initC]*1000)
    initC['natal_kick_1'] = np.random.uniform(0, 100, 1000)
    initC['phi_1'] = np.random.uniform(-90, 90, 1000)
    initC['theta_1'] = np.random.uniform(0, 360, 1000)
    initC['mean_anomaly_1'] = np.random.uniform(0, 360, 1000)
    initC['porb'] = np.random.uniform(50, 190, 1000)
    initC['sep'] = utils.a_from_p(p=initC.porb.values, m1=initC.mass_1.values, m2=initC.mass_2.values)
    initC['bin_num'] = np.linspace(0, 1000, 1000)

    bpp_restart, bcm_restart, initC_restart, kick_info_restart = Evolve.evolve(
        initialbinarytable=initC, BSEDict=BSEDict, SSEDict=SSEDict
    )

    bpp_BH = bpp_restart.loc[
        (bpp_restart.kstar_1 == 14)
        & (bpp_restart.kstar_2 == 1)
        & (bpp_restart.porb > 0)].groupby('bin_num', as_index=False).first()

    bpp_BH[['tphys', 'mass_1', 'mass_2', 'porb', 'ecc']]
