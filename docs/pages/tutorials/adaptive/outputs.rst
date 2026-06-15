.. _adaptive_outputs:

***************************************
Handling outputs from adaptive sampling
***************************************

This tutorial assumes that you've already gone through :ref:`adaptive_basics`.

Reading your results from a file
================================

After you've finished running your adaptive sampling simulation, you will now have some results stored as :class:`~cosmic.output.COSMICStroopOutput` object. If you saved these results to a file, then you can reload them by running

.. code-block:: python

    from cosmic.output import COSMICStroopOutput

    results = COSMICStroopOutput.from_file("path/to/your/file.h5")

Understanding your outputs
==========================

The :class:`~cosmic.output.COSMICStroopOutput` class stores all of the information you need to analyse your simulation. Let's step through some of the different attributes that you will need, and assume for the purposes of this guide that you have ``N`` samples, in ``D`` dimensions, with ``H`` hits.

Sample information
------------------

Each :class:`~cosmic.output.COSMICStroopOutput` object contains a full record of every sample that was made during the simulation. 

- ``samples``: contains a every sampled point and as such has shape ``(N, D)``
- ``param_names``: is a list of length ``D`` with names corresponding to each column in the ``samples`` array
- ``is_hit``: is an array of length ``N`` with boolean values for whether a sample was a hit (and as such will sum to ``H``)

Hit details
-----------

For the actual hits (i.e. the samples that you most care about), this class also stores the full evolution history. In particular, ``bpp``, ``bcm``, ``initC``, and ``kick_info`` all contain the usual ``COSMIC`` evolution tables (see :ref:`evolve_single` if you're not familiar).

.. tip::

    If you want to just explore your hits in particular, you can sub-select them by doing

    .. code-block:: python

        just_hits = results[results.is_hit]

    which masks the class just like you would with a :class:`~cosmic.output.COSMICOutput`. Be aware you likely still need the full population for access to the weights (we'll cover weights below).

General metadata
----------------

The class also stores other metadata that you may find useful. These include:

- ``num_explored``, which is the number of systems evolved during the exploration phase
- ``num_hits``, which is the total raw hit count across all phases
- ``fraction_explored``, which is the fraction of total systems used for exploration.


How to interpret adaptive sampling weights
==========================================

So now for the important intuition part. ``COSMIC`` and ``STROOPWAFEL`` have now provided you with a sample of a rare population by preferentially sampling the parameter space that you've specified. However, we of course want to account for the fact that this *is* a rare population. This is where the weights come in. Each sample is assigned an adaptive importance sampling weight, which tells you how rare this sample is (smaller weights are rarer).

.. admonition:: TODO

    Add maths

This means that if you want to plot a true distribution of your sampled systems -- let's say the primary mass -- you need to use the weights in your plotting.

.. code-block:: python

    import matplotlib.pyplot as plt
    from cosmic.output import COSMICStroopOutput

    results = COSMICStroopOutput.from_file("YOUR_SIMULATION.h5")
    primary_mass = results.samples[:, results.param_names.index("mass_1")]

    plt.hist(primary_mass, weights=results.weights, bins=50, density=True)
    plt.xlabel(r"Primary mass, $m_1$ [$\rm M_\odot$]")
    plt.ylabel("Probability density")
    plt.show()


.. warning::

    You should **always** apply your weights to any plot that you make. In histograms you can supply them directly. For scatter plots you could consider changing the size of your points, or using a 2D histogram or ``hexbin`` instead.

Drawing a representative sample from your simulation
====================================================