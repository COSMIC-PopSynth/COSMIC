*****************************************
Evolving a Monte-Carlo sampled population
*****************************************
Once an initial binary population is sampled, it can be evolved using the ``Evolve`` class just as we've done so far.
You can read more about sampling initial binary populations in the :ref:`runpop` page.

Note that the same process used other examples applies here as well: the ``BSEDict`` must be supplied,
if the flags in the dictionary change from their defaults.

First, let's import the necessary modules from COSMIC:

.. ipython:: python

    from cosmic.sample.initialbinarytable import InitialBinaryTable

    from cosmic.evolve import Evolve


Now let's use the independent sampler to generate an initial binary population. You can learn more 
about the independent sampler in the :ref:`independent` page.

.. ipython:: python

    final_kstars = range(16)

    InitialBinaries, mass_singles, mass_binaries, n_singles, n_binaries = \
          InitialBinaryTable.sampler('independent', final_kstars, final_kstars,
                                     binfrac_model=0.5, primary_model='kroupa01',
                                     ecc_model='sana12', porb_model='sana12',
                                     qmin=-1, SF_start=13700.0,
                                     SF_duration=0.0, met=0.02, size=10, keep_singles=True)


And finally, we can evolve the initial binary population using the Evolve class as we've done in the previous
guides:

.. include:: ../../../_generated/default_bsedict.rst

.. ipython:: python
    :okwarning:

    bpp, bcm, initC, kick_info = Evolve.evolve(
        initialbinarytable=InitialBinaries, BSEDict=BSEDict
    )

    print(bcm.iloc[:10])

    print(bpp)
