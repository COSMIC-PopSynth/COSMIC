.. _independent:

*************************
Independent distributions
*************************

This guide will show you how to sample an initial binary population using the independent sampler in COSMIC.
First import the :class:`~cosmic.sample.initialbinarytable.InitialBinaryTable` class and the independent sampler

.. ipython::

    In [1]: from cosmic.sample.initialbinarytable import InitialBinaryTable

    In [2]: from cosmic.sample.sampler import independent

The independent sampler contains multiple models for each binary parameter.
You can access the available models using the independent sampler help call:

.. ipython::

    In [3]: help(independent.get_independent_sampler)


Targetting specific final kstar types
=====================================

The ``final_kstar1`` and ``final_kstar2`` parameters are lists that contain the kstar types
that you would like the final population to contain.

The final kstar is the final state of the binary system we are interested in and is based on the BSE kstar naming convention, see :ref:`kstar-table` for more information.

Thus, if you want to generate a
population containing double white dwarfs with CO and ONe WD primaries and He-WD secondaries,
the final kstar inputs would be:

.. ipython::

    In [4]: final_kstar1 = [11, 12]

    In [5]: final_kstar2 = [10]

Since we are interested in binaries, we only retain the binary systems that are likely to produce the user specified final kstar types.
However, we also keep track of the total mass of the single and binary stars as well as the number of binary and single stars so that we can scale our results to larger populations.
If you don't want to filter the binaries, you can supply final kstars as

.. ipython::

    In [6]: final_kstars = np.linspace(0, 14, 15)

    In [7]: InitialBinaries, mass_singles, mass_binaries, n_singles, n_binaries = InitialBinaryTable.sampler('independent', final_kstars, final_kstars, binfrac_model=0.5, primary_model='kroupa01', ecc_model='sana12', porb_model='sana12', qmin=-1, SF_start=13700.0, SF_duration=0.0, met=0.02, size=10000)

Additionally if you are interested in single stars then you can specify ``keep_singles=True``. In this case, the singles will be added onto the end of the InitialBinaryTable where ``kstar_1`` will host the singles, ``kstar_2`` will be filled with 15s only, and all orbital properties (e.g. ``porb`` or ``ecc``) will be indicated with -1.

Understanding parameter sampling models
=======================================

Similar to the help for the sampler, the different models that can be used for each parameter
to be sampled can be accessed by the help function for the argument. The syntax for each parameter
sample is always: sample_`parameter`. See the example for the star formation
history (SFH) below:

.. ipython::

    In [8]: help(independent.Sample.sample_SFH)

Stopping conditions for sampling
================================

In addition to telling ``COSMIC`` *how* to sample, you also need to tell it how *much* to sample. You can
do this in one of two ways: (1) by specifying the number of binaries to sample, or (2) by specifying the total mass of singles and binaries to sample. Let's look at both of these in more detail.

Number of binaries
------------------

Using the final kstar inputs we mentioned above, the initial binary population can be sampled based on the desired number of binaries (10000 in this case) as follows

.. ipython::

    In [9]: InitialBinaries, mass_singles, mass_binaries, n_singles, n_binaries = InitialBinaryTable.sampler('independent', final_kstar1, final_kstar2, binfrac_model=0.5, primary_model='kroupa01', ecc_model='sana12', porb_model='sana12', qmin=-1, SF_start=13700.0, SF_duration=0.0, met=0.02, size=10000)

    In [10]: print(InitialBinaries)

.. note::
    
    The length of the initial binary data set, ``InitialBinaries``, does not always match
    the size parameter provided to :meth:`~cosmic.sample.initialbinarytable.InitialBinaryTable.sampler`.
    This is because of the various cuts that the sampler makes to the population (e.g. the binary fraction,
    which is either a fraction between 0 and 1 or mass dependent following the
    prescription in `van Haaften+2013 <http://adsabs.harvard.edu/abs/2012A%26A...537A.104V>`_.) specified by the user.

Total mass sampled
------------------

Alternatively, we could do the same thing but now instead set our ``sampling_target`` to be the total mass and aim for 15000 solar masses. This is done by setting ``sampling_target="total_mass"`` and ``total_mass=15000``.

.. ipython::

    In [10]: InitialBinaries, mass_singles, mass_binaries, n_singles, n_binaries = InitialBinaryTable.sampler('independent', final_kstar1, final_kstar2, binfrac_model=0.5, primary_model='kroupa01', ecc_model='sana12', porb_model='sana12', qmin=-1, SF_start=13700.0, SF_duration=0.0, met=0.02, sampling_target="total_mass", total_mass=15000)

    In [11]: print(InitialBinaries)

And we can check what the total sampled mass was by looking at the sum of the ``mass_singles`` and ``mass_binaries`` variables

.. ipython::

    In [12]: print(mass_singles + mass_binaries)

.. tip::

    If you'd like to avoid your sample overshooting your desired ``total_mass`` and instead get as close to this value as possible,
    you can set ``trim_extra_samples=True``. This will trim the sample to get a total mass as close as possible to your target.
    In many cases, this will be within a solar mass, but could be as large as twice the maximum stellar mass (for the very rare case that
    the final binary drawn is the most massive primary with an equal mass ratio).

Mass dependent binary fractions and mass pairings
=================================================

If you want to have separate binary fractions and mass pairings for low and high mass stars, you can by supplying the ``msort`` kwarg to the sampler. This sets the mass above which an alternative mass pairing (specified by kwargs ``qmin_msort`` and ``m2_min_msort``) and binary fraction model (specified by kwarg ``binfrac_model_msort``) are used. This is handy if you want, for example, a higher binary fraction and more equal mass pairings for high mass stars.

Below we show the effect of different assumptions for the independent initial sampler. The standard assumptions are shown in purple, the assumptions of `Sana et al. 2012 <https://ui.adsabs.harvard.edu/abs/2012Sci...337..444S/abstract>`_ are shown in orange, and the assumptions of `Moe et al. 2019 <https://ui.adsabs.harvard.edu/abs/2019ApJ...875...61M/abstract>`_ are shown in green.

.. plot::
   :include-source: False

    >>> from cosmic.utils import a_from_p
    >>> from cosmic.sample.initialbinarytable import InitialBinaryTable
    >>> import pandas as pd
    >>> import numpy as np
    >>> import matplotlib.pyplot as plt
    >>> n_samples = 10000
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
    >>>                                                                                         size=n_samples)
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
    >>>                                                                                         size=n_samples)
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
    >>>                                                                                    size=n_samples)
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
