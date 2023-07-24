.. _runpop:

####################################
Generate a binary population by hand
####################################

The process to generate a synthetic binary population, is similar to the
process to evolve a single/multiple binaries by hand: first generate an
initial population, then evolve it with the Evolve class.

An initialized binary population consists of a collection of binary systems
with assigned primary and secondary masses, orbital periods, eccentricities,
metallicities, and star formation histories. These parameters are randomly
sampled from observationally motivated distribution functions.

In COSMIC, the initial sample is done through an initial binary sampler which works
with the InitialBinaryTable class. There are two samplers available:

1. `independent` : initialize binaries with independent parameter
distributions for the primary mass, mass ratio, eccentricity, separation,
and binary fraction

2. `multidim` : initialize binaries with multidimensional parameter
distributions according to `Moe & Di Stefano 2017 <http://adsabs.harvard.edu/abs/2017ApJS..230...15M>`_

We consider both cases below.

***********
independent
***********

First import the InitialBinaryTable class and the independent sampler

.. ipython::

    In [1]: from cosmic.sample.initialbinarytable import InitialBinaryTable

    In [2]: from cosmic.sample.sampler import independent

The independent sampler contains multiple models for each binary parameter.
You can access the available models using the independent sampler help call:

.. ipython::

    In [3]: help(independent.get_independent_sampler)

The final_kstar1 and final_kstar2 parameters are lists that contain the kstar types
that you would like the final population to contain.

The final kstar is the final state of the binary system we are interested in and is based on the BSE kstar naming convention, see :ref:`kstar-table` for more information.

Thus, if you want to generate a
population containing double white dwarfs with CO and ONe WD primaries and He-WD secondaries,
the final kstar inputs would be:

.. ipython::

    In [4]: final_kstar1 = [11, 12]

    In [5]: final_kstar2 = [10]

Similar to the help for the sampler, the different models that can be used for each parameter
to be sampled can be accessed by the help function for the argument. The syntax for each parameter
sample is always: sample_`parameter`. See the example for the star formation
history (SFH) below:

.. ipython::

    In [6]: help(independent.Sample.sample_SFH)

Using the final kstar inputs above, the initial binary population is sampled as:

.. ipython::

    In [6]: InitialBinaries, mass_singles, mass_binaries, n_singles, n_binaries = InitialBinaryTable.sampler('independent', final_kstar1, final_kstar2, binfrac_model=0.5, primary_model='kroupa01', ecc_model='sana12', porb_model='sana12', qmin=-1, m2_min=0.08, SF_start=13700.0, SF_duration=0.0, met=0.02, size=10000)

    In [7]: print(InitialBinaries)

NOTE: the length of the initial binary data set, InitialBinaries, does not always match
the size parameter provided to InitialBinaryTable.sampler.
This is because the sampler accounts for a binary fraction specified by the user with the binfrac_model parameter, which is either a fraction between 0 and 1 or mass dependent following the prescription in `van Haaften+2013 <http://adsabs.harvard.edu/abs/2012A%26A...537A.104V>`_.

If you want to have separate binary fractions and mass pairings for low and high mass stars, you can by supplying the `msort` kwarg to the sampler. This sets the mass above which an alternative mass pairing (specified by kwargs `qmin_msort` and `m2_min_msort`) and binary fraction model (specified by kwarg `binfrac_model_msort`) are used. This is handy if you want, for example, a higher binary fraction and more equal mass pairings for high mass stars.

Since we are interested in binaries, we only retain the binary systems that are likely to produce the user specified final kstar types. However, we also keep track of the total mass of the single and binary stars as well as the number of binary and single stars so that we can scale our results to larger populations. If you don't want to filter the binaries, you can supply final kstars as

.. ipython::

    In [8]: final_kstars = np.linspace(0, 14, 15)

    In [9]: InitialBinaries, mass_singles, mass_binaries, n_singles, n_binaries = InitialBinaryTable.sampler('independent', final_kstars, final_kstars, binfrac_model=0.5, primary_model='kroupa01', ecc_model='sana12', porb_model='sana12', qmin=-1, m2_min=0.08, SF_start=13700.0, SF_duration=0.0, met=0.02, size=10000)

Below we show the effect of different assumptions for the independent initial sampler. The standard assumptions are shown in blue, the assumptions of `Sana et al. 2012 <https://ui.adsabs.harvard.edu/abs/2012Sci...337..444S/abstract>`_ are shown in orange, and the assumptions of `Moe et al. 2019 <https://ui.adsabs.harvard.edu/abs/2019ApJ...875...61M/abstract>`_.

.. plot::
   :include-source: False

    >>> from cosmic.utils import a_from_p
    >>> from cosmic.sample.initialbinarytable import InitialBinaryTable
    >>> import pandas as pd
    >>> import numpy as np
    >>> import matplotlib.pyplot as plt
    >>> final_kstar = np.linspace(0,14,15)
    >>> colors = {'green' : '#1b9e77', 'orange' : '#d95f02', 'purple' : '#7570b3'}
    >>> initC_logP, m_sin_logP, m_bin_logP, n_sin_logP, n_bin_logP = InitialBinaryTable.sampler('independent',
    >>>                                                                                         final_kstar1=final_kstar,
    >>>                                                                                         final_kstar2=final_kstar,
    >>>                                                                                         binfrac_model=1.0,
    >>>                                                                                         primary_model='kroupa01',
    >>>                                                                                         ecc_model='thermal',
    >>>                                                                                         porb_model='log_uniform',
    >>>                                                                                         qmin=-1,
    >>>                                                                                         SF_start=13700.0,
    >>>                                                                                         SF_duration=0.0,
    >>>                                                                                         met=0.02,
    >>>                                                                                         size=100000)
    >>> initC_Sana, m_sin_Sana, m_bin_Sana, n_sin_Sana, n_bin_Sana = InitialBinaryTable.sampler('independent',
    >>>                                                                                         final_kstar1=final_kstar,
    >>>                                                                                         final_kstar2=final_kstar,
    >>>                                                                                         binfrac_model=1.0,
    >>>                                                                                         primary_model='kroupa01',
    >>>                                                                                         ecc_model='sana12',
    >>>                                                                                         porb_model='sana12',
    >>>                                                                                         qmin=-1,
    >>>                                                                                         SF_start=13700.0,
    >>>                                                                                         SF_duration=0.0,
    >>>                                                                                         met=0.02,
    >>>                                                                                         size=100000)
    >>> initC_Moe, m_sin_Moe, m_bin_Moe, n_sin_Moe, n_bin_Moe = InitialBinaryTable.sampler('independent',
    >>>                                                                                    final_kstar1=final_kstar,
    >>>                                                                                    final_kstar2=final_kstar,
    >>>                                                                                    binfrac_model=1.0,
    >>>                                                                                    primary_model='kroupa01',
    >>>                                                                                    ecc_model='sana12',
    >>>                                                                                    porb_model='moe19',
    >>>                                                                                    qmin=-1,
    >>>                                                                                    SF_start=13700.0,
    >>>                                                                                    SF_duration=0.0,
    >>>                                                                                    met=0.02,
    >>>                                                                                    size=100000)
    >>> 
    >>> initC_logP['sep'] = a_from_p(p=initC_logP.porb, m1=initC_logP.mass_1, m2=initC_logP.mass_2)
    >>> initC_Sana['sep'] = a_from_p(p=initC_Sana.porb, m1=initC_Sana.mass_1, m2=initC_Sana.mass_2)
    >>> initC_Moe['sep'] = a_from_p(p=initC_Moe.porb, m1=initC_Moe.mass_1, m2=initC_Moe.mass_2)
    >>> fig = plt.figure(figsize = (15,6))
    >>> ax1 = plt.subplot(231)
    >>> ax2 = plt.subplot(232)
    >>> ax3 = plt.subplot(233)
    >>> ax4 = plt.subplot(234)
    >>> ax5 = plt.subplot(235)
    >>> ax6 = plt.subplot(236)
    >>> ax1.hist(np.log10(initC_logP.mass_1), bins = 20, histtype='step', density=True,
    >>>          lw=3, color=colors['purple'], label='independent')
    >>> ax1.hist(np.log10(initC_Sana.mass_1), bins = 20, histtype='step', density=True,
    >>>          lw=3, color=colors['orange'], label='Sana+2012')
    >>> ax1.hist(np.log10(initC_Moe.mass_1), bins = 20, histtype='step', density=True,
    >>>          lw=3, color=colors['green'], label='Moe+2019')
    >>> ax1.set_xlabel(r'Log$_{10}$(M$_1$/M$_{\odot}$)', size=18)
    >>> ax1.set_ylabel('normalized counts', size=18)
    >>> ax1.legend(prop={'size' : 18})
    >>> ax2.hist(np.log10(initC_logP.porb), bins = 20, histtype='step', density=True,
    >>>          lw=3, color=colors['purple'], label='independent')
    >>> ax2.hist(np.log10(initC_Sana.porb), bins = 20, histtype='step', density=True,
    >>>          lw=3, color=colors['orange'], label='Sana+2012')
    >>> ax2.hist(np.log10(initC_Moe.porb), bins = 20, histtype='step', density=True,
    >>>          lw=3, color=colors['green'], label='Moe+2019')
    >>> ax2.set_xlabel(r'Log$_{10}$(P$_{\rm{orb}}$/day)', size=18)
    >>> ax3.hist(initC_logP.ecc, bins = 10, histtype='step', density=True,
    >>>          lw=3, color=colors['purple'], label='independent')
    >>> ax3.hist(initC_Sana.ecc, bins = 10, histtype='step', density=True,
    >>>          lw=3, color=colors['orange'], label='Sana+2012')
    >>> ax3.hist(initC_Moe.ecc, bins = 10, histtype='step', density=True,
    >>>          lw=3, color=colors['green'], label='Moe+2019')
    >>> ax3.set_xlabel('Eccentricity', size=18)
    >>> ax4.hist(initC_logP.mass_2/initC_logP.mass_1, bins = 20, histtype='step', density=True,
    >>>          lw=3, color=colors['purple'], label='independent')
    >>> ax4.hist(initC_Sana.mass_2/initC_Sana.mass_1, bins = 20, histtype='step', density=True,
    >>>          lw=3, color=colors['orange'], label='Sana+2012')
    >>> ax4.hist(initC_Moe.mass_2/initC_Moe.mass_1, bins = 20, histtype='step', density=True,
    >>>          lw=3, color=colors['green'], label='Moe+2019')
    >>> ax4.set_xlabel(r'q=M$_1$/M$_2$', size=18)
    >>> ax4.set_ylabel('normalized counts', size=18)
    >>> ax5.hist(np.log10(initC_logP.sep), bins = 20, histtype='step', density=True,
    >>>          lw=3, color=colors['purple'], label='independent')
    >>> ax5.hist(np.log10(initC_Sana.sep), bins = 20, histtype='step', density=True,
    >>>          lw=3, color=colors['orange'], label='Sana+2012')
    >>> ax5.hist(np.log10(initC_Moe.sep), bins = 20, histtype='step', density=True,
    >>>          lw=3, color=colors['green'], label='Moe+2019')
    >>> ax5.set_xlabel(r'Log$_{10}$(a/R$_{\odot}$)', size=18)
    >>> ax6.hist(np.log10(initC_logP.sep*(1-initC_logP.ecc)), bins = 20, histtype='step', density=True,
    >>>          lw=3, color=colors['purple'], label='independent')
    >>> ax6.hist(np.log10(initC_Sana.sep*(1-initC_Sana.ecc)), bins = 20, histtype='step', density=True,
    >>>          lw=3, color=colors['orange'], label='Sana+2012')
    >>> ax6.hist(np.log10(initC_Moe.sep*(1-initC_Moe.ecc)), bins = 20, histtype='step', density=True,
    >>>          lw=3, color=colors['green'], label='Moe+2019')
    >>> ax6.set_xlabel(r'Log$_{10}$(a(1-e)/R$_{\odot}$)', size=18)
    >>> fig.tight_layout()
    >>> fig.show()


********
multidim
********

COSMIC implements multidimensionally distributed initial binaries according to `Moe & Di Stefano 2017 <http://adsabs.harvard.edu/abs/2017ApJS..230...15M>`_. The python code used in COSMIC to create this sample was written by Mads Sorenson, and is based on the IDL codes written to accompany `Moe & Di Stefano 2017 <http://adsabs.harvard.edu/abs/2017ApJS..230...15M>`_.

The multidimensional initial binary data is sampled in COSMIC as follows:

.. ipython::

    In [1]: from cosmic.sample.initialbinarytable import InitialBinaryTable

    In [2]: from cosmic.sample.sampler import multidim

To see the arguments necessary to call the multidimensional sampler use the help function:

.. ipython::

    In [3]: help(multidim.get_multidim_sampler)

The random seed is used to reproduce your initial sample, since there are several stochastic processes involved in the muldimensional sample.
As in the independent sampler, the final_kstar1 and final_kstar2 inputs are lists containing the kstar types that the evolved population should contain.

The multidimensional sample is generated as follows:

.. ipython::

    In [4]: InitialBinaries, mass_singles, mass_binaries, n_singles, n_binaries = InitialBinaryTable.sampler('multidim', final_kstar1=[11], final_kstar2=[11], rand_seed=2, nproc=1, SF_start=13700.0, SF_duration=0.0, met=0.02, size=10)

    In [5]: print(InitialBinaries)

.. note::

    NOTE that in the multidimensional case, the binary fraction is a parameter in the sample. This results in the size of the initial binary data matching the size provided to the sampler. As in the independent sampling case, we keep track of the total sampled mass of singles and binaries as well as the total number of single and binary stars to scale thesimulated population to astrophysical populations.

.. plot::
   :include-source: False

    >>> from cosmic.utils import a_from_p
    >>> from cosmic.sample.initialbinarytable import InitialBinaryTable
    >>> import pandas as pd
    >>> import numpy as np
    >>> import matplotlib.pyplot as plt
    >>> final_kstar = np.linspace(0,14,15)
    >>> colors = {'green' : '#1b9e77', 'purple' : '#d95f02', 'orange' : '#7570b3'}
    >>> final_kstar = np.linspace(0,14,15)
    >>> initC_mult, m_sin_mult, m_bin_mult, n_sin_mult, n_bin_mult = InitialBinaryTable.sampler('multidim',
    >>>                                                                                         final_kstar1=final_kstar,
    >>>                                                                                         final_kstar2=final_kstar,
    >>>                                                                                         rand_seed=2,
    >>>                                                                                         nproc=1,
    >>>                                                                                         SF_start=13700.0,
    >>>                                                                                         SF_duration=0.0,
    >>>                                                                                         met=0.02,
    >>>                                                                                         size=100000)
    >>> initC_mult['sep'] = a_from_p(p=initC_mult.porb, m1=initC_mult.mass_1, m2=initC_mult.mass_2)

*************************************
Evolving an initial binary population
*************************************
As in :ref:`examples`, once an initial binary population is sampled, it is evolved using the Evolve class. Note that the same process used in :ref:`examples` applies here as well: the BSEDict must be supplied, but only need be resupplied if the flags in the dictionary change.

The syntax for the Evolve class is as follows:

.. ipython::
    :okwarning:

    In [1]: from cosmic.evolve import Evolve

    In [2]: BSEDict = {'xi': 1.0, 'bhflag': 1, 'neta': 0.5, 'windflag': 3, 'wdflag': 1, 'alpha1': 1.0, 'pts1': 0.001, 'pts3': 0.02, 'pts2': 0.01, 'epsnov': 0.001, 'hewind': 0.5, 'ck': 1000, 'bwind': 0.0, 'lambdaf': 0.0, 'mxns': 3.0, 'beta': -1.0, 'tflag': 1, 'acc2': 1.5, 'grflag' : 1, 'remnantflag': 4, 'ceflag': 0, 'eddfac': 1.0, 'ifflag': 0, 'bconst': 3000, 'sigma': 265.0, 'gamma': -2.0, 'pisn': 45.0, 'natal_kick_array' : [[-100.0,-100.0,-100.0,-100.0,0.0], [-100.0,-100.0,-100.0,-100.0,0.0]], 'bhsigmafrac' : 1.0, 'polar_kick_angle' : 90, 'qcrit_array' : [0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0], 'cekickflag' : 2, 'cehestarflag' : 0, 'cemergeflag' : 0, 'ecsn' : 2.25, 'ecsn_mlow' : 1.6, 'aic' : 1, 'ussn' : 0, 'sigmadiv' :-20.0, 'qcflag' : 1, 'eddlimflag' : 0, 'fprimc_array' : [2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0], 'bhspinflag' : 0, 'bhspinmag' : 0.0, 'rejuv_fac' : 1.0, 'rejuvflag' : 0, 'htpmb' : 1, 'ST_cr' : 1, 'ST_tide' : 1, 'bdecayfac' : 1, 'rembar_massloss' : 0.5, 'kickflag' : 0, 'zsun' : 0.019, 'bhms_coll_flag' : 0, 'don_lim' : -1, 'acc_lim' : -1, 'rtmsflag' : 0}

    In [3]: bpp, bcm, initC, kick_info  = Evolve.evolve(initialbinarytable=InitialBinaries, BSEDict=BSEDict)

    In [4]: print(bcm.iloc[:10])

    In [5]: print(bpp)


*****************
ClusterMonteCarlo
*****************

New in COSMIC 3.4, you can now use COSMIC to sample initial conditions that can be used in the simulation of a Globular Cluster (GC), using the ClusterMonteCarlo (CMC) software package. To create these initial conditions, and save them in a format readable by CMC, you can do the following.

.. ipython::

    In [1]: from cosmic.sample.initialcmctable import InitialCMCTable

    In [2]: from cosmic.sample.sampler import cmc

To see the arguments necessary to call the CMC sampler use the help function:

.. ipython::

    In [3]: help(cmc.get_cmc_sampler)

.. ipython::
    :okwarning:

    In [1]: from cosmic.sample import InitialCMCTable
      
    In [2]: Singles, Binaries = InitialCMCTable.sampler('cmc', binfrac_model=0.2, primary_model='kroupa01', ecc_model='sana12', porb_model='sana12', qmin=-1.0, cluster_profile='plummer', met=0.014, size=40000, params='../examples/Params.ini', gamma=4, r_max=100)

    In [3]: InitialCMCTable.write(Singles, Binaries, filename="input.hdf5")

    In [4]: InitialCMCTable.write(Singles, Binaries, filename="input.fits")
