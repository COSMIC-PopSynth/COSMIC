.. _adaptive_basics:

*********************************************
Adaptive importance sampling for rare systems
*********************************************

The standard :ref:`independent <independent>` and :ref:`multidimensional <multidim>`
samplers draw binary initial parameters from distributions that are defined for each sampler. These samples can then be evolved with COSMIC to find the final population. For many scenarios this works well: common-envelope episodes, mass-transferring
binaries, and white dwarf systems all occur frequently enough that thousands of random draws
yield a workable sample.

However, some outcomes are extremely rare. Binary black holes that merge within the Hubble time
form at very low rates (even more so for NS + NS mergers, especially at high metallicity).
Sampling such populations with regular Monte Carlo draws may require tens of millions of COSMIC calls
to accumulate even a few hundred systems — which may end up being prohibitive for a large parameter survey.

``COSMIC`` includes a vectorised implementation of the ``STROOPWAFEL`` algorithm
(`Broekgaarden et al. 2019 <https://doi.org/10.1093/mnras/stz2558>`_) that solves this
problem using *adaptive importance sampling*.  The sampler first explores parameter space
to locate the progenitor regions of the target population, then concentrates its simulation
budget on those regions, and finally corrects the biased sampling with importance weights so
that any weighted statistic remains an unbiased estimator of the true prior-weighted
distribution.

This implementation is based on the implementation designed by Lokesh Khandelwal, Floris Kummer, and Stephen Justham, which built upon the original STROOPWAFEL algorithm.


When should I use this?
=======================

Use ``STROOPWAFEL`` whenever you need a statistically representative sample of a rare binary
outcome and cannot afford the total binary count that a flat Monte Carlo draw would require.
Example use cases include:

* Double black holes (or neutron stars) merging within the Hubble time (GW sources)
* Long lived BH + stellar-companion systems

If your target population is common the plain independent sampler is simpler and fast enough.

How it works: think Battleships
===============================

``STROOPWAFEL`` works in three phases that you can think of as a game of Battleships. You wouldn't continue to shoot randomly at the grid after you find a ship — instead, you would concentrate your fire around the hit location to sink it. Similarly, ``STROOPWAFEL`` first explores parameter space with random draws, then adapts a proposal distribution based on the hits it finds, and finally concentrates its sampling from that proposal.

**Exploration — random fire**
    Binaries are drawn at random from the prior distributions and evolved with ``COSMIC``.
    Every binary that produces the desired outcome is recorded as a *hit*. An adaptive
    stopping criterion (based on the observed hit rate) decides when the remaining budget
    is better spent on refinement than on further random exploration.

**Adaptation — mark the ships**
    One multivariate Gaussian component is placed at each hit location in parameter space.
    The width of each Gaussian is derived from the local density of the prior (via the CDF
    of the prior distribution) so that the proposal is appropriately broad regardless of
    the parameter's scale. Together the Gaussians form a *mixture model* — a coarse map of
    where progenitors live.

**Refinement — concentrate fire**
    New binaries are drawn from the Gaussian mixture instead of the broad prior.  Because
    the mixture is concentrated near known progenitor regions, a much larger fraction of the
    simulated systems produce hits.  Between refinement generations the mixture is
    optionally updated via an expectation-maximisation (EM) step that shifts weight toward
    the most productive Gaussians.

**Weight calculation — correct the bias**
    Every simulated system receives an importance weight

    .. math::

        w(x) = \frac{\pi(x)}{Q(x)}, \qquad
        Q(x) = f_e\,\pi(x) + (1 - f_e)\,q(x),

    where :math:`\pi(x)` is the prior probability density, :math:`q(x)` is the Gaussian
    mixture density, and :math:`f_e` is the fraction of systems drawn from the prior
    during exploration. Using these weights, any statistic computed on the hit population
    is an unbiased estimator of the corresponding prior-weighted quantity.

Setup
=====

Define the parameter space
--------------------------

The :class:`~cosmic.sample.stroopwafel.ParameterSpace` class defines which binary
parameters are sampled and their distributions.  Each
:class:`~cosmic.sample.stroopwafel.Parameter` specifies a name, physical-space bounds, and
a distribution.  The distribution sets both how the parameter is drawn and its prior
probability density — in importance sampling these are one and the same.

.. code-block:: python

    import numpy as np
    from cosmic.sample.stroopwafel import ParameterSpace, Parameter

    params = ParameterSpace([
        Parameter('mass_1',      5.0,        150.0,      dist='kroupa'),
        Parameter('q',           0.01,       1.0,        dist='uniform'),
        Parameter('porb',        10**(0.15), 10**(5.5),  dist='sana'),
        Parameter('ecc',         1e-9,       0.9999,     dist='sana_ecc'),
        Parameter('metallicity', 0.0001,     0.03,       dist='flat_in_log'),
    ])

The ``dist`` argument names one of the built-in distributions. Several common initial-distribution choices are available by default - the Kroupa IMF, Sana orbital periods and eccentricities, flat-in-log metallicity, and more — and lets you define your own just as easily.  See :ref:`adaptive_distributions` for the full list and how to add custom distributions.

Parameters are stored and returned in the order you provide them, which sets the column order of every sample array. Use ``params.names`` to check the order and ``params.idx('param_name')`` to get the index of a particular parameter.

.. note::

    Bounds are always given in **physical** space.  The ``'sana'`` period is sampled in
    :math:`\log_{10}(P / \mathrm{day})`, so the example writes its bounds as ``10**(0.15)``
    and ``10**(5.5)`` - the distribution applies the :math:`\log_{10}` transform internally.


Complete the binary definition
------------------------------

Every binary handed to ``COSMIC`` is defined by **five** parameters: ``mass_1``,
``mass_2``, ``porb``, ``ecc``, and ``metallicity``.

Each one must be provided in **exactly one** of two ways: either it is sampled (a
:class:`~cosmic.sample.stroopwafel.Parameter` with that name in your ``ParameterSpace``) or
it is returned by an optional ``derive_params`` function.  If any of the five is neither
sampled nor derived, :class:`~cosmic.sample.stroopwafel.AdaptiveSampler` raises an error as soon as it is constructed.

In the parameter space above we sampled ``mass_1``, ``porb``, ``ecc``, and ``metallicity``
directly, but sampled the mass *ratio* ``q`` rather than ``mass_2``.  ``derive_params``
fills in the gap.  It receives a dictionary mapping each sampled name to its ``(N,)`` array
of physical values and returns a dictionary of the remaining parameters (a scalar is
broadcast to all ``N`` binaries):

.. code-block:: python

    def derive_params(sampled):
        """Provide binary parameters not drawn from the ParameterSpace."""
        return {'mass_2': sampled['mass_1'] * sampled['q']}

This design makes it easy to adaptively sample in just a few dimensions while holding the
rest fixed.  For instance, to explore *only* orbital period you would sample ``porb`` and
fix everything else:

.. code-block:: python

    params = ParameterSpace([
        Parameter('porb', 10**(0.15), 10**(5.5), dist='sana'),
    ])

    def derive_params(sampled):
        return {'mass_1': 30.0, 'mass_2': 25.0, 'ecc': 0.0, 'metallicity': 0.02}

If all five parameters are sampled directly, ``derive_params`` is unnecessary and can be
omitted entirely.


Choose a rejection function
---------------------------

Before a batch is passed to ``COSMIC``, unphysical systems are filtered out: stars already
overflowing their Roche lobes at ZAMS, binaries whose components are in contact, and
systems below the hydrogen-burning limit for the secondary.  The built-in
:func:`~cosmic.sample.stroopwafel.rejection.default_reject` function performs all of these
checks:

It receives the assembled binary parameters as a dictionary (``mass_1``, ``mass_2``,
``porb``, ``ecc``, ``metallicity``) and returns a boolean mask where ``True`` means the
system is rejected; the orbital separation is computed internally from ``porb``.  For most
use cases this default is appropriate.  If you need extra cuts — for example imposing a
minimum secondary mass — you can wrap it:

.. code-block:: python

    def my_reject(binary_params):
        base_mask = default_reject(binary_params)
        # additionally reject secondaries below 1 M_sun
        base_mask |= (binary_params['mass_2'] < 1.0)
        return base_mask

Passing ``reject_systems`` is optional; you can instead pass ``None`` to skip physical rejection entirely.


Identify what constitutes a hit
-------------------------------

The ``is_interesting`` argument identifies which evolved systems count as hits.  It receives
the ``COSMIC`` ``bpp`` DataFrame for the current batch and must return a tuple
``(n_hits, hit_bin_nums)`` where ``hit_bin_nums`` is an integer array of ``bin_num`` values
(0-indexed within the batch).

The ``STROOPWAFEL`` sampler comes with two preset functions in :mod:`cosmic.sample.stroopwafel.presets`.
These functions focus on double compact objects (DCOs) and are suitable for gravitational-wave source studies:

``any_dco(kstar_1, kstar_2)``
    Selects all bound double compact objects whose stellar types match the supplied lists,
    regardless of merger time.

``merging_dco(kstar_1, kstar_2, max_merge_time=13.7)``
    Like ``any_dco`` but additionally requires the merger time to be
    less than ``max_merge_time`` Gyr.

.. note::

    The ``merging_dco`` function requires the `LEGWORK python package <https://legwork.readthedocs.io/en/latest/>`_ to compute the merger time.  If you do not have LEGWORK installed, ``merging_dco`` will raise an error.

.. code-block:: python

    from cosmic.sample.stroopwafel.presets import any_dco, merging_dco

    # All bound BH-BH systems, no merger time cut
    is_interesting = any_dco(kstar_1=[14], kstar_2=[14])

    # Only NS-NS systems merging within the Hubble time
    is_interesting_merging = merging_dco(kstar_1=[13], kstar_2=[13], max_merge_time=13.7)


You can also write a fully custom hit function.  For example, to find BH + stellar companion systems that remain bound for at least 100 Myr after the BH forms:

.. code-block:: python

    import numpy as np

    def bh_star_100myr(bpp):
        """Hit: BH with a stellar companion bound for at least 100 Myr."""
        bh_star = bpp.loc[
            (
                ((bpp['kstar_1'] == 14) & (bpp['kstar_2'] < 10))
                | ((bpp['kstar_2'] == 14) & (bpp['kstar_1'] < 10))
            )
            & (bpp['sep'] > 0)
        ]
        if bh_star.empty:
            return 0, np.array([], dtype=int)

        # Duration in the BH + star state for each binary
        span = bh_star.groupby('bin_num')['tphys'].agg(lambda t: t.max() - t.min())
        hits = span.index[span >= 100.0].values
        return len(hits), hits


Running the sampler
===================

Now we can put it all together! You can run the sampler with the main :class:`~cosmic.sample.stroopwafel.AdaptiveSampler` class.  The most important arguments are the parameter space, the total number of systems to evolve, the batch size, the BSE and SSE physics settings, the hit function, the ``derive_params`` function (if needed), and the rejection function. See the API documentation (:class:`~cosmic.sample.stroopwafel.AdaptiveSampler`) for a full list of options.

Let's try this out with a few examples.

Examples
--------

Bound BH + BH binaries
^^^^^^^^^^^^^^^^^^^^^^

Let's imagine we want to sample the population of bound BH + BH binaries. We can use the same parameter space and ``derive_params`` function as above, and the built-in ``any_dco`` hit function to select all bound BH + BH systems. We can sample over a five-dimensional parameter space covering primary mass, mass ratio, orbital period, eccentricity, and metallicity.

First we can import the necessary parts from the ``cosmic.sample.stroopwafel`` module, the preset hit function, and setup a BSEDict.

.. code-block:: python

    import numpy as np
    from cosmic.sample.stroopwafel import AdaptiveSampler, ParameterSpace, Parameter
    from cosmic.sample.stroopwafel.presets import any_dco

.. include:: ../../../_generated/default_bsedict.rst

COSMIC v4+ also expects an ``SSEDict`` of single stellar evolution settings (which selects
the stellar engine).  We'll use the default ``sse`` engine here, you can swap in METISSE too if you like!

.. code-block:: python

    SSEDict = {'stellar_engine': 'sse'}

Then we can define a simple parameter space, where we avoid sampling low-mass primaries since we know they
cannot produce a BH.

.. code-block:: python

    params = ParameterSpace([
        Parameter('mass_1',      5.0,        150.0,      dist='kroupa'),
        Parameter('q',           0.01,       1.0,        dist='uniform'),
        Parameter('porb',        10**(0.15), 10**(5.5),  dist='sana'),
        Parameter('ecc',         1e-9,       0.9999,     dist='sana_ecc'),
        Parameter('metallicity', 0.0001,     0.03,       dist='flat_in_log'),
    ])

Since we only sampled the mass ratio ``q``, we need to derive the secondary mass from the primary mass and ``q``:

.. code-block:: python

    def derive_params(sampled):
        return {'mass_2': sampled['mass_1'] * sampled['q']}


And then it's just a matter of setting it going!

.. code-block:: python

    sampler = AdaptiveSampler(
        parameter_space=params,
        total_systems=50_000,           # adjust this for more samples
        batch_size=500,                 # adjust this to sample more or fewer systems per call to COSMIC
        BSEDict=BSEDict,
        SSEDict=SSEDict,
        is_interesting=any_dco(kstar_1=[14], kstar_2=[14]),
        derive_params=derive_params,
        reject_systems="default",
        nproc=4,
        n_generations=1,
        seed=42,
    )
    result = sampler.run()

    print(f"Total hits:        {result.num_hits}")
    print(f"Weighted hit rate: {result.hit_rate:.4e} ± {result.hit_rate_uncertainty:.4e}")


BH + star binaries surviving 100 Myr
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Now let's repeat that whole scenario, but instead of BH + BH binaries we want to sample BH + stellar companion systems that remain bound for at least 100 Myr after the BH forms. We can use the same parameter space and ``derive_params`` function as above, but this time we will use the custom ``bh_star_100myr`` hit function defined earlier.

.. code-block:: python

    import numpy as np
    from cosmic.sample.stroopwafel import AdaptiveSampler
    from cosmic.sample.stroopwafel.rejection import default_reject

    sampler = AdaptiveSampler(
        parameter_space=params,          # reuse from BHBH example
        total_systems=50_000,
        batch_size=500,
        BSEDict=BSEDict,                 # reuse from BHBH example
        SSEDict=SSEDict,                 # reuse from BHBH example
        is_interesting=bh_star_100myr,   # we defined this earlier
        derive_params=derive_params,     # reuse from BHBH example
        reject_systems="default",
        nproc=4,
        n_generations=1,
        seed=42,
    )
    result = sampler.run()


Rules of thumb
==============

Choosing ``total_systems``, ``batch_size``, and ``n_generations`` is something of an art but there are some rules of thumb to get you started.  The optimal settings depend on the rarity and complexity of the target population, the dimensionality of the parameter space, and your computational resources.

``batch_size``
--------------

``batch_size`` sets how many systems are passed to :meth:`~cosmic.evolve.Evolve.evolve` per call.

* Aim for ``batch_size`` to be a multiple of ``nproc`` so that COSMIC distributes work
  evenly across cores.
* Values of 200-1000 are typical.  Batches smaller than ~50 increase Python overhead per
  call; batches larger than ~5000 may cause memory pressure from the output DataFrames.

``total_systems``
-----------------

This is the total number of binary evolutions across all phases. It's hard to know how many you'll need without knowing the rarity of the target population. As a rule of thumb, aim for at least ~30 hits during exploration before the adaptation phase begins — fewer hits lead to a poorly-constrained Gaussian mixture.  If exploration ends with very few hits, increase ``total_systems`` and re-run.

``n_generations``
-----------------

Each refinement generation uses an equal share of the remaining budget after exploration.
The EM step between generations can improve the mixture, but with diminishing returns. You are probably safe with just 1 generation unless you have a very rare population.

.. tip::

    To run a plain Monte Carlo without any adaptation or refinement — useful as a baseline
    or for common populations — pass ``mc_only=True``.  The sampler will draw from the
    prior only, and importance weights will equal the prior density divided by itself
    (i.e. all weights are equal).

Saving your results
===================

Once you have your samples, you can save them to disk as an HDF5 file with the :meth:`~cosmic.output.COSMICStroopOutput.save` method.  This saves the parameter samples, derived quantities, and hit information that can be loaded later for analysis.

.. code-block:: python

    result.save('bhbh_samples.h5')

We'll talk more about how to load and analyse these results in the :ref:`adaptive_outputs` tutorial! But first, let's look at how to define custom distributions for your parameters in the next tutorial :ref:`adaptive_distributions`.