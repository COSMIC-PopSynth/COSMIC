.. _adaptive_outputs:

***************************************
Handling outputs from adaptive sampling
***************************************

This tutorial assumes that you've already gone through :ref:`adaptive_basics`.

In this tutorial we're going to cover how to read in and interpret the outputs from an adaptive sampling run.
We'll also cover how to draw a representative sample from your simulation, and how to use the weights that are generated during the adaptive sampling process.

Reading your results from a file
================================

After you've finished running your adaptive sampling simulation, you will now have some results stored as :class:`~cosmic.output.COSMICStroopOutput` object.
If you saved these results to a file, then you can reload them by running

.. code-block:: python

    from cosmic.output import COSMICStroopOutput

    results = COSMICStroopOutput.from_file("path/to/your/file.h5")

Understanding your outputs
==========================

The :class:`~cosmic.output.COSMICStroopOutput` class stores all of the information you need to analyse your simulation.
Let's step through some of the different attributes that you will need, and assume for the purposes of this guide that you have ``N`` samples, in ``D`` dimensions, with ``H`` hits.

Sample information
------------------

Each :class:`~cosmic.output.COSMICStroopOutput` object contains a full record of every sample that was made during the simulation. 

- ``samples``: contains a every sampled point and as such has shape ``(N, D)``
- ``param_names``: is a list of length ``D`` with names corresponding to each column in the ``samples`` array
- ``is_hit``: is an array of length ``N`` with boolean values for whether a sample was a hit (and as such will sum to ``H``)

Hit details
-----------

For the actual hits (i.e. the samples that you most care about), this class also stores the full evolution history.
In particular, ``bpp``, ``bcm``, ``initC``, and ``kick_info`` all contain the usual ``COSMIC`` evolution tables (see :ref:`evolve_single` if you're not familiar).

.. tip::

    If you want to just explore your hits in particular, you can sub-select them by doing

    .. code-block:: python

        just_hits = results[results.is_hit]

    which masks the class just like you would with a :class:`~cosmic.output.COSMICOutput` (see :ref:`analysis_interface` if you're not familiar).
    Be aware you likely still need the full population for access to the weights (we'll cover weights below).

General metadata
----------------

The class also stores other metadata that you may find useful. These include:

- ``num_explored``, which is the number of systems evolved during the exploration phase
- ``num_hits``, which is the total raw hit count across all phases
- ``fraction_explored``, which is the fraction of total systems used for exploration.

In addition, two derived properties summarise the weighted population:

- ``hit_rate``, the importance-weighted hit rate (an unbiased estimate of the true rate)
- ``hit_rate_uncertainty``, the standard error on that rate


How to interpret adaptive sampling weights
==========================================

So now for the important intuition part. ``COSMIC`` and ``STROOPWAFEL`` have now provided you with a sample of a rare population by preferentially sampling the parameter space that you've specified. However, we of course want to account for the fact that this *is* a rare population. This is where the weights come in. Each sample is assigned an adaptive importance sampling weight, which tells you how rare this sample is (smaller weights are rarer).

Concretely, each system :math:`x` is assigned an importance weight

.. math::

    w(x) = \frac{\pi(x)}{Q(x)}, \qquad
    Q(x) = f_e\,\pi(x) + (1 - f_e)\,q(x),

where :math:`\pi(x)` is the prior (astrophysical) probability density, :math:`q(x)` is the
density of the adapted Gaussian mixture, and :math:`f_e` is the fraction of the budget spent
on exploration.  The denominator :math:`Q(x)` is therefore the *actual* density from which
the system was drawn — a blend of the broad prior (during exploration) and the concentrated
mixture (during refinement).

In words, the weight is the ratio of how often a system *should* appear under the prior to
how often it *actually* appeared under the combined sampling scheme.  A hit discovered deep
inside an oversampled region picks up a small weight (we drew far more of them than nature
would), while a system drawn straight from the prior has a weight close to one.  Summed over
the population these weights turn any statistic into an unbiased estimator of the
corresponding prior-weighted quantity; for example ``sum(weights[is_hit]) / N`` is an
unbiased estimate of the true hit rate, which is exactly what
:attr:`~cosmic.output.COSMICStroopOutput.hit_rate` returns.

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

Applying weights at plot time is the right approach for visualising distributions, but
sometimes you want an actual *set of systems* that is representative of the underlying
population. For example, you may want a fixed number of binaries for further analysis that don't require any weights.
Because the simulation deliberately oversamples the rare region, you cannot take the hits at
face value; you need to resample them in proportion to their weights.

We provide a convenience method for this in :class:`~cosmic.output.COSMICStroopOutput`, but the underlying procedure is simple. It is a standard weighted bootstrap: draw indices from the hit population with probability proportional to their weights, with replacement.

.. code-block:: python

    import numpy as np
    from cosmic.output import COSMICStroopOutput

    results = COSMICStroopOutput.from_file("YOUR_SIMULATION.h5")
    representative_sample, bin_nums = results.draw_representative_sample(n_samples=1000)

The resulting ``representative_sample`` array provides a set of 1000 systems with their parameters drawn from the underlying population, and ``bin_nums`` provides the corresponding indices into the original population so that you can access the full evolution history if you need it.

Wrap-up
=======

And that's everything you need to know about working with the outputs from your adaptive sampling run. You can now read in your results, understand the weights, and draw a representative sample from your simulation.

Finally, we're going to look at how to checkpoint your adaptive sampling run so that you can resume it later in :ref:`adaptive_checkpoint` - see you there!
