*********************************
Evolving a fixed grid of binaries
*********************************

Sometimes it is helpful to run a fixed grid of initial binaries to explore how
changing a single parameter affects the evolved binary. Let's start by importing the necessary modules:

.. ipython::

    In [1]: from cosmic.sample.initialbinarytable import InitialBinaryTable

    In [2]: from cosmic.evolve import Evolve


Here we evolve the same system that produces a GW150914-like binary, but run over several initial orbital
periods spaced evenly in log space.

.. ipython::
    :okwarning:

    In [3]: n_grid = 10

    In [4]: binary_grid = InitialBinaryTable.InitialBinaries(m1=np.ones(n_grid)*100.0,
       ...:                                                  m2=np.ones(n_grid)*85.0,
       ...:                                                  porb=np.logspace(3,5,n_grid),
       ...:                                                  ecc=np.ones(n_grid)*0.65,
       ...:                                                  tphysf=np.ones(n_grid)*13700.0,
       ...:                                                  kstar1=np.ones(n_grid),
       ...:                                                  kstar2=np.ones(n_grid),
       ...:                                                  metallicity=np.ones(n_grid)*0.005)

    In [5]: BSEDict = {'xi': 1.0, 'bhflag': 1, 'neta': 0.5, 'windflag': 3, 'wdflag': 1, 'alpha1': 1.0,
       ...:            'pts1': 0.001, 'pts3': 0.02, 'pts2': 0.01, 'epsnov': 0.001, 'hewind': 0.5, 'ck': 1000,
       ...:            'bwind': 0.0, 'lambdaf': 0.0, 'mxns': 3.0, 'beta': -1.0, 'tflag': 1, 'acc2': 1.5,
       ...:            'grflag' : 1, 'remnantflag': 4, 'ceflag': 0, 'eddfac': 1.0, 'ifflag': 0,
       ...:            'bconst': 3000, 'sigma': 265.0, 'gamma': -2.0, 'pisn': 45.0,
       ...:            'natal_kick_array' : [[-100.0,-100.0,-100.0,-100.0,0.0], [-100.0,-100.0,-100.0,-100.0,0.0]],
       ...:            'bhsigmafrac' : 1.0, 'polar_kick_angle' : 90,
       ...:            'qcrit_array' : [0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],
       ...:            'cekickflag' : 2, 'cehestarflag' : 0, 'cemergeflag' : 0, 'ecsn' : 2.25,
       ...:            'ecsn_mlow' : 1.6, 'aic' : 1, 'ussn' : 0, 'sigmadiv' :-20.0, 'qcflag' : 1, 'eddlimflag' : 0,
       ...:            'fprimc_array' : [2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0],
       ...:            'bhspinflag' : 0, 'bhspinmag' : 0.0, 'rejuv_fac' : 1.0, 'rejuvflag' : 0, 'htpmb' : 1,
       ...:            'ST_cr' : 1, 'ST_tide' : 1, 'bdecayfac' : 1, 'rembar_massloss' : 0.5, 'kickflag' : 1,
       ...:            'zsun' : 0.014, 'bhms_coll_flag' : 0, 'don_lim' : -1, 'acc_lim' : -1, 'rtmsflag' : 0,
       ...:            'wd_mass_lim': 1}

    In [6]: print(binary_grid)

    In [7]: bpp, bcm, initC, kick_info = Evolve.evolve(initialbinarytable=binary_grid, BSEDict=BSEDict)

    In [8]: print(bpp)

    In [9]: print(bcm)