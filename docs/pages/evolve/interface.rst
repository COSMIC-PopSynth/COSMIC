**************************
Analysing your simulations
**************************

After running your COSMIC simulations with an ``Evolve.evolve()`` call, you will have several outputs to analyse.
The outputs are:

- ``bpp`` : A pandas DataFrame containing the binary population properties at important timesteps for each binary.
- ``bcm`` : A pandas DataFrame containing information about binaries at user-defined timesteps
- ``initC`` : A pandas DataFrame containing the initial conditions of each binary.
- ``kick_info`` : A pandas DataFrame containing information about natal kicks imparted to compact objects during supernovae.

These outputs can be combined into a single ``COSMICOutput`` object for easier analysis and plotting.

Creating a ``COSMICOutput`` object
==================================

You can create a ``COSMICOutput`` object by passing in the outputs from your evolution. Let's start by 
samples ~100 binaries and evolving them.

.. ipython:: python

    from cosmic.sample.initialbinarytable import InitialBinaryTable
    from cosmic.evolve import Evolve
    from cosmic.output import COSMICOutput
    import matplotlib.pyplot as plt
    import numpy as np

    InitialBinaries, mass_singles, mass_binaries, n_singles, n_binaries = InitialBinaryTable.sampler(
        'independent', [13, 14], [13, 14], binfrac_model=0.5, primary_model='kroupa01',
        ecc_model='sana12', porb_model='sana12', qmin=-1, SF_start=13700.0, SF_duration=0.0,
        met=0.002, size=1000)

    BSEDict = {
        'xi': 1.0, 'bhflag': 1, 'neta': 0.5, 'windflag': 3, 'wdflag': 1, 'alpha1': 1.0,
        'pts1': 0.001, 'pts3': 0.02, 'pts2': 0.01, 'epsnov': 0.001, 'hewind': 0.5, 'ck': 1000,
        'bwind': 0.0, 'lambdaf': 0.0, 'mxns': 3.0, 'beta': -1.0, 'tflag': 1, 'acc2': 1.5,
        'grflag' : 1, 'remnantflag': 4, 'ceflag': 0, 'eddfac': 1.0, 'ifflag': 0, 'bconst': 3000,
        'sigma': 265.0, 'gamma': -2.0, 'pisn': 45.0,
        'natal_kick_array' : [[-100.0,-100.0,-100.0,-100.0,0.0], [-100.0,-100.0,-100.0,-100.0,0.0]],
        'bhsigmafrac' : 1.0, 'polar_kick_angle' : 90,
        'qcrit_array' : [0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],
        'cekickflag' : 2, 'cehestarflag' : 0, 'cemergeflag' : 0, 'ecsn' : 2.25, 'ecsn_mlow' : 1.6,
        'aic' : 1, 'ussn' : 0, 'sigmadiv' :-20.0, 'qcflag' : 1, 'eddlimflag' : 0,
        'fprimc_array' : [2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0],
        'bhspinflag' : 0, 'bhspinmag' : 0.0, 'rejuv_fac' : 1.0, 'rejuvflag' : 0, 'htpmb' : 1,
        'ST_cr' : 1, 'ST_tide' : 1, 'bdecayfac' : 1, 'rembar_massloss' : 0.5, 'kickflag' : 5,
        'zsun' : 0.014, 'bhms_coll_flag' : 0, 'don_lim' : -1, 'acc_lim' : -1, 'rtmsflag' : 0,
        'wd_mass_lim': 1
    }

    bpp, bcm, initC, kick_info = Evolve.evolve(initialbinarytable=InitialBinaries, BSEDict=BSEDict)

Now we can create a ``COSMICOutput`` object quite easily (with an optional label):

.. ipython:: python

    output = COSMICOutput(bpp=bpp, bcm=bcm, initC=initC, kick_info=kick_info,
                          label="My First COSMIC Output")
    print(output)

This object now contains all the relevant data from our simulation in one place, and we can use its methods to analyse and plot the results.


Saving and loading from a file
==============================

You can save a ``COSMICOutput`` object to a file for later use or sharing, and load it back when needed. The
entire file gets saved to a single HDF5 file.

.. ipython:: python

    output.save('your_output_file.h5')

    print(output.initC.head())

    # Later, or in another script
    loaded_output = COSMICOutput(file='your_output_file.h5')

    print(loaded_output.initC.head())


.. note::

    This file will also save the version of COSMIC used to create the output. The load function will then
    warn you if you are loading an output created with a different version of COSMIC than the one you are
    currently using.


Plotting population distributions
=================================

The ``COSMICOutput`` class includes several built-in plotting functions to visualize the results of your simulations.
First, let's plot the initial mass distribution of the primary stars in our binaries.

.. ipython:: python
    :okwarning:
    :okexcept:
    
    fig, ax = output.plot_distribution(x_col="mass_1", when="initial", show=False);
    @savefig initial_mass_distribution.png
    plt.show()

.. image:: initial_mass_distribution.png
    :alt:
    :width: 100%

We could also compare this to the final mass distribution of the primary stars after evolution.

.. ipython:: python
    :okwarning:
    :okexcept:

    output.plot_distribution(x_col="mass_1", when="final", show=False);
    @savefig final_mass_distribution.png
    plt.show()

.. image:: final_mass_distribution.png
    :alt:
    :width: 100%

In addition to histograms, you can create scatter plots to visualize relationships between different parameters.

.. ipython:: python
    :okwarning:
    :okexcept:

    fig, ax = output.plot_distribution(
        x_col="teff_1", y_col="lum_1", c_col="mass_1",
        when="final", show=False
    );
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.invert_xaxis()
    
    @savefig hrd.png
    plt.show()

.. image:: hrd.png
    :alt:
    :width: 100%

And we also can colour the points by the stellar type of the primary star at the end of evolution and get a custom
colour map for these ones.

.. ipython:: python
    :okwarning:
    :okexcept:

    fig, ax = output.plot_distribution(
        x_col="mass_1", y_col="porb", c_col="kstar_1",
        when="final", show=False
    );
    ax.set_xscale("log")
    ax.set_yscale("log")
    
    @savefig m1_porb_kstar.png
    plt.show()

.. image:: m1_porb_kstar.png
    :alt:
    :width: 100%

Any column name from the ``bpp`` (or ``initC`` for initial conditions) DataFrames can be used for the x, y, and colour col names.

Subselecting binaries
=====================

Often you'll want to focus on a specific subset of binaries from your simulation. 
You can easily subselect binaries simply by indexing the ``COSMICOutput`` object. You can index based on the
specific ``bin_num`` of the binaries, or by using a boolean mask.

Let's say you just want the first binary, or the binaries with ``bin_num`` 0, 2, and 4, or even every binary
between 10 and 20, you just do:

.. ipython:: python

    first_binary = output[0]
    some_binaries = output[[0, 2, 4]]
    range_of_binaries = output[10:21]

But perhaps you only care about binaries where the primary star starts with at least 5 solar masses and the
binary merges during evolution. You can create a boolean mask and use that to index the ``COSMICOutput`` object:

.. ipython:: python

    mask = (output.initC['mass_1'] >= 5.0) & (output.final_bpp["sep"] == 0.0)
    selected_binaries = output[mask]
    print(selected_binaries)


Re-running with more detailed output
====================================

Now let's say you care about one binary in particular, and you want to re-run the evolution of just that binary with more detailed output.
You can do this easily by subselecting the binary you care about, and then calling the ``rerun_with_settings()``
method on the resulting ``COSMICOutput`` object with a smaller ``dtp``.

Let's find a binary that forms a black hole and re-run it with more detailed output:

.. ipython:: python

    bh = output[(output.final_bpp["kstar_1"] == 14)
                | (output.final_bpp["kstar_2"] == 14)]
    detailed_bh = bh.rerun_with_settings(new_settings={'dtp': 0.0}, inplace=False)
    print(detailed_bh.bcm)

Now clearly staring at the DataFrame isn't very helpful, so let's plot the evolution of the binary's separation over time.

.. ipython:: python
    :okwarning:
    :okexcept:

    bh_bin_num = detailed_bh.initC['bin_num'].iloc[0]
    detailed_bh.plot_detailed_evolution(bin_num=bh_bin_num, show=False);
    @savefig detailed_bh.png
    plt.show()

.. image:: detailed_bh.png
    :alt:
    :width: 100%

This shows the full evolution, but we can ignore any time a while after the binary forms a BH
by setting a maximum time for the x-axis.

.. ipython:: python
    :okwarning:
    :okexcept:

    formed_a_bh = (detailed_bh.bcm['kstar_1'] == 14) | (detailed_bh.bcm['kstar_2'] == 14)
    t_max = detailed_bh.bcm['tphys'][formed_a_bh].min() + 10.0  # 10 Myr after first BH forms

    detailed_bh.plot_detailed_evolution(bin_num=bh_bin_num,
                                        t_max=t_max, show=False);
    @savefig detailed_bh_limited.png
    plt.show()

.. image:: detailed_bh_limited.png
    :alt:
    :width: 100%


Re-running with new physics settings
====================================

You can also re-run the evolution of your binaries with different physics settings, but the same initial conditions.
This is useful for testing how different assumptions impact your results. You can do this using the ``rerun_with_settings()`` method, passing in a dictionary of
new settings. For example, let's say we want to see how changing the common envelope efficiency affects our results.

.. ipython:: python

    ce_alpha_10 = output.rerun_with_settings(
        new_settings={'alpha1': 10}, inplace=False
    )
    n_merger_original = len(output.final_bpp[output.final_bpp["sep"] == 0.0])
    n_merger_ce_alpha_10 = len(ce_alpha_10.final_bpp[ce_alpha_10.final_bpp["sep"] == 0.0])
    print(f"Original number of mergers: {n_merger_original}")
    print(f"Number of mergers with ceflag=10: {n_merger_ce_alpha_10}")

This shows how making common-envelope evolution more efficient allows more binaries to survive and avoid merging.


Additionally, if your new settings affect natal kicks (e.g., changing remnant mass prescriptions), you can 
also choose to reset the natal kicks when re-running the evolution. This ensures that the kicks are sampled
according to the new physics settings. You can do this by setting the ``reset_kicks`` parameter to ``True`` in the ``rerun_with_settings()`` method.

.. ipython:: python

    kickflag_1 = output.rerun_with_settings(
        new_settings={'kickflag': 1}, reset_kicks=True, inplace=False
    )
    print(f"Average kick velocity with original settings: {output.kick_info['natal_kick'].mean():1.2f} km/s")
    print(f"Average kick velocity with kickflag=1 and reset kicks: {kickflag_1.kick_info['natal_kick'].mean():1.2f} km/s")


Examining output from `cosmic-pop`
==================================

If you have run a population synthesis simulation using ``cosmic-pop``, you can also load the output
from that simulation into a ``COSMICPopOutput`` object. This object extends the functionality of ``COSMICOutput``
to handle the binary and single star populations from ``cosmic-pop``, as well as storing information about the
convergence of the simulation.

Let's say that you've saved a ``cosmic-pop`` simulation to a file called ``dat_kstar1_13_14_kstar2_13_14_SFstart_13700.0_SFduration_0.0_metallicity_0.02.h5``.
You can load this file into a ``COSMICPopOutput`` object like so:

.. code-block:: python

    from cosmic.output import COSMICPopOutput

    pop_output = COSMICPopOutput(file='dat_kstar1_13_14_kstar2_13_14_SFstart_13700.0_SFduration_0.0_metallicity_0.02.h5')
    print(pop_output)
    

This file contains both binary and single star populations, which you can access via the ``output`` and ``singles_output`` attributes, respectively.
These are both ``COSMICOutput`` objects, so you can use all the same methods and attributes on them as well (e.g., plotting distributions, subselecting binaries, re-running with new settings, etc.).

.. code-block:: python

    # the full bpp table for binaries
    print(pop_output.output.bpp)

If you set ``keep_singles=False`` when running ``cosmic-pop``, the ``singles_output`` attribute will be ``None``.

You can access the usual cosmic-pop outputs as attributes of the class, such as:

.. code-block:: python

    print(pop_output.conv)
    print(pop_output.match)
    print(pop_output.n_binaries)
    print(pop_output.mass_stars)

Combining binary and single star outputs
----------------------------------------

If you want to combine the binary and single star outputs into a single ``COSMICOutput`` object, you can use the ``to_combined_output()`` method.
This will concatenate the binary and single star DataFrames together. This can let you more easily analyse the full stellar population from your simulation.

.. code-block:: python

    combined_output = pop_output.to_combined_output()
    
    # for example, we could select the BHBH binaries and re-run them with more detailed output
    bhbh_binaries = combined_output[
        (combined_output.final_bpp["kstar_1"] == 14) &
        (combined_output.final_bpp["kstar_2"] == 14)
    ]
    detailed_bhbh = bhbh_binaries.rerun_with_settings(new_settings={'dtp': 0.1}, inplace=False)
    detailed_bhbh.plot_detailed_evolution(
        bin_num=detailed_bhbh.initC['bin_num'].iloc[0],
        t_max=100
    );


Wrapping up
===========

The ``COSMICOutput`` class provides a convenient way to manage and analyse the results of your COSMIC simulations.
With built-in methods for plotting, subselecting binaries, and re-running evolutions with different settings,
it will hopefully make your analysis workflow smoother and get you to your interesting results faster!