******************************
Multidimensional distributions
******************************

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
    >>>                                                                                         size=10000)
    >>> initC_mult['sep'] = a_from_p(p=initC_mult.porb, m1=initC_mult.mass_1, m2=initC_mult.mass_2)

