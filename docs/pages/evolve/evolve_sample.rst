*****************************************
Evolving a Monte-Carlo sampled population
*****************************************
Once an initial binary population is sampled, it can be evolved using the ``Evolve`` class just as we've done so far.
You can read more about sampling initial binary populations in the :ref:`runpop` page.

Note that the same process used other examples applies here as well: the ``BSEDict`` must be supplied,
if the flags in the dictionary change from their defaults.

First, let's import the necessary modules from COSMIC:

.. ipython::

    In [1]: from cosmic.sample.initialbinarytable import InitialBinaryTable

    In [2]: from cosmic.evolve import Evolve


Now let's use the independent sampler to generate an initial binary population. You can learn more 
about the independent sampler in the :ref:`independent` page.

.. ipython::

    In [3]: final_kstars = range(16)

    In [4]: InitialBinaries, mass_singles, mass_binaries, n_singles, n_binaries = \
       ...:       InitialBinaryTable.sampler('independent', final_kstars, final_kstars,
       ...:                                  binfrac_model=0.5, primary_model='kroupa01',
       ...:                                  ecc_model='sana12', porb_model='sana12',
       ...:                                  qmin=-1, SF_start=13700.0,
       ...:                                  SF_duration=0.0, met=0.02, size=10, keep_singles=True)


And finally, we can evolve the initial binary population using the Evolve class as we've done in the previous
guides:

.. ipython::
    :okwarning:

    In [5]: BSEDict = {'xi': 1.0, 'bhflag': 1, 'neta': 0.5, 'windflag': 3, 'wdflag': 1, 'alpha1': 1.0, 'pts1': 0.001, 'pts3': 0.02, 'pts2': 0.01, 'epsnov': 0.001, 'hewind': 0.5, 'ck': 1000, 'bwind': 0.0, 'lambdaf': 0.0, 'mxns': 3.0, 'beta': -1.0, 'tflag': 1, 'acc2': 1.5, 'grflag' : 1, 'remnantflag': 4, 'ceflag': 0, 'eddfac': 1.0, 'ifflag': 0, 'bconst': 3000, 'sigma': 265.0, 'gamma': -2.0, 'pisn': 45.0, 'natal_kick_array' : [[-100.0,-100.0,-100.0,-100.0,0.0], [-100.0,-100.0,-100.0,-100.0,0.0]], 'bhsigmafrac' : 1.0, 'polar_kick_angle' : 90, 'qcrit_array' : [0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0], 'cekickflag' : 2, 'cehestarflag' : 0, 'cemergeflag' : 0, 'ecsn' : 2.25, 'ecsn_mlow' : 1.6, 'aic' : 1, 'ussn' : 0, 'sigmadiv' :-20.0, 'qcflag' : 1, 'eddlimflag' : 0, 'fprimc_array' : [2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0], 'bhspinflag' : 0, 'bhspinmag' : 0.0, 'rejuv_fac' : 1.0, 'rejuvflag' : 0, 'htpmb' : 1, 'ST_cr' : 1, 'ST_tide' : 1, 'bdecayfac' : 1, 'rembar_massloss' : 0.5, 'kickflag' : 1, 'zsun' : 0.019, 'bhms_coll_flag' : 0, 'don_lim' : -1, 'acc_lim' : -1, 'rtmsflag' : 0, 'wd_mass_lim' : 1}

    In [6]: bpp, bcm, initC, kick_info = Evolve.evolve(initialbinarytable=InitialBinaries,
       ...:                                            BSEDict=BSEDict)

    In [7]: print(bcm.iloc[:10])

    In [8]: print(bpp)