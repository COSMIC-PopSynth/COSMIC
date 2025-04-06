*******************
Restarting a binary
*******************

COSMIC allows you to restart a binary from any point in its evolution from a COSMIC generated bpp array.

.. ipython::

    In [1]: from cosmic.sample.initialbinarytable import InitialBinaryTable

    In [2]: from cosmic.evolve import Evolve

Below we provide an example of the same evolutionary track
started from the beginning and three different points in the evolution:

- sometime between the beginning and the first object going supernova
- between the first and second supernova
- after both supernova

.. ipython::
    :okwarning:

    In [3]: single_binary = InitialBinaryTable.InitialBinaries(m1=25.543645, m2=20.99784, porb=446.795757, ecc=0.448872, tphysf=13700.0, kstar1=1, kstar2=1, metallicity=0.002)
    
    In [4]: BSEDict = {'xi': 1.0, 'bhflag': 1, 'neta': 0.5, 'windflag': 3, 'wdflag': 1, 'alpha1': 5.0, 'pts1': 0.001, 'pts3': 0.02, 'pts2': 0.01, 'epsnov': 0.001, 'hewind': 0.5, 'ck': 1000, 'bwind': 0.0, 'lambdaf': 0.0, 'mxns': 3.0, 'beta': -1.0, 'tflag': 1, 'acc2': 1.5, 'remnantflag': 3, 'ceflag': 0, 'eddfac': 1.0, 'ifflag': 0, 'bconst': 3000, 'sigma': 265.0, 'gamma': -2.0, 'pisn': 45.0, 'natal_kick_array' : [[-100.0,-100.0,-100.0,-100.0,0.0], [-100.0,-100.0,-100.0,-100.0,0.0]], 'bhsigmafrac' : 1.0, 'polar_kick_angle' : 90, 'qcrit_array' : [0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0], 'cekickflag' : 2, 'cehestarflag' : 0, 'cemergeflag' : 0, 'ecsn' : 2.5, 'ecsn_mlow' : 1.4, 'aic' : 1, 'ussn' : 0, 'sigmadiv' :-20.0, 'qcflag' : 1, 'eddlimflag' : 0, 'fprimc_array' : [2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0], 'bhspinflag' : 0, 'bhspinmag' : 0.0, 'rejuv_fac' : 1.0, 'rejuvflag' : 0, 'htpmb' : 1, 'ST_cr' : 1, 'ST_tide' : 0, 'bdecayfac' : 1, 'randomseed' : -1235453, 'grflag' : 1, 'rembar_massloss' : 0.5, 'kickflag' : 1, 'zsun' : 0.014,  'grflag' : 1, 'bhms_coll_flag' : 0, 'don_lim' : -1, 'acc_lim' : -1, 'rtmsflag' : 0, 'wd_mass_lim': 1}
    
    In [5]: for i in [3, 7, 11]:
       ...:     bpp, bcm, initC, kick_info = Evolve.evolve(initialbinarytable=single_binary, BSEDict=BSEDict)
       ...:     for column in bpp.columns:
       ...:         initC = initC.assign(**{column:bpp.iloc[i][column]})
       ...:     bpp_mid, bcm_mid, initC_mid, kick_info = Evolve.evolve(initialbinarytable=initC, BSEDict={})
       ...:     if i == 3:
       ...:         print("From beginning")
       ...:         print(bpp)
       ...:     print("Started in middle at Index {0}".format(i))
       ...:     print(bpp_mid)


Natal kick example
==================

One example of where restarting a binary can be extremely helpful is studying
how natal kicks affect a binary independently of its previous evolution. This is
particularly relevant for
`Gaia BH1 <https://ui.adsabs.harvard.edu/abs/2023MNRAS.518.1057E/abstract>`_
and `Gaia BH2 <https://ui.adsabs.harvard.edu/abs/2023MNRAS.521.4323E/abstract>`_ which are difficult to produce through the standard common envelope
channels. We can still study the effect of natal kicks on these binaries if we
restart the evolution after the mass transfer would occur. We can do this by using a binary which gets us to the right masses given the metallicity, then overwrite some of the initial conditions to resample the natal kicks and pre-explosion separation. 

.. ipython::
    :okwarning:

    In [6]: from cosmic import utils
       ...: import pandas as pd

    In [7]: single_binary = InitialBinaryTable.InitialBinaries(m1=65.0, m2=0.93, porb=4500, ecc=0.448872, 
       ...:                                                    tphysf=13700.0, kstar1=1, kstar2=1, metallicity=0.014*0.6)

    In [8]: bpp, bcm, initC, kick_info = Evolve.evolve(initialbinarytable=single_binary, BSEDict=BSEDict)

    In [9]: for column in bpp.columns:
       ...:     initC = initC.assign(**{column:bpp.iloc[6][column]})

    In [10]: initC = pd.concat([initC]*1000)
       ....: initC['natal_kick_1'] = np.random.uniform(0, 100, 1000)
       ....: initC['phi_1'] = np.random.uniform(-90, 90, 1000)
       ....: initC['theta_1'] = np.random.uniform(0, 360, 1000)
       ....: initC['mean_anomaly_1'] = np.random.uniform(0, 360, 1000)
       ....: initC['porb'] = np.random.uniform(50, 190, 1000)
       ....: initC['sep'] = utils.a_from_p(p=initC.porb.values, m1=initC.mass_1.values, m2=initC.mass_2.values)
       ....: initC['bin_num'] = np.linspace(0, 1000, 1000)

    In [11]: bpp_restart, bcm_restart, initC_restart, kick_info_restart = Evolve.evolve(initialbinarytable=initC, BSEDict={})

    In [12]: bpp_BH = bpp_restart.loc[(bpp_restart.kstar_1 == 14) & (bpp_restart.kstar_2 == 1) & (bpp_restart.porb > 0)].groupby('bin_num', as_index=False).first()

    In [13]: bpp_BH[['tphys', 'mass_1', 'mass_2', 'porb', 'ecc']]

