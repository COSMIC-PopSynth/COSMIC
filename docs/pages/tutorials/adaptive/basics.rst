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


When should I use this?
=======================

Use ``STROOPWAFEL`` whenever you need a statistically representative sample of a rare binary
outcome and cannot afford the total binary count that flat Monte Carlo would require.
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
:class:`~cosmic.sample.stroopwafel.Parameter` specifies a name, physical-space bounds, a
sampling distribution, and a prior distribution.

.. code-block:: python

    import numpy as np
    from cosmic.sample.stroopwafel import ParameterSpace, Parameter

    params = ParameterSpace([
        Parameter('mass_1',      5.0,    150.0,  sampler='kroupa',     prior='kroupa'),
        Parameter('q',           0.01,   1.0,    sampler='uniform',    prior='uniform'),
        Parameter('porb',        0.15,   5.5,    sampler='sana',       prior='sana'),
        Parameter('ecc',         1e-9,   0.9999, sampler='sana_ecc',   prior='sana_ecc'),
        Parameter('metallicity', 0.0001, 0.03,   sampler='flat_in_log',prior='flat_in_log'),
    ])

Parameters are stored and returned in alphabetical order by name. Use ``params.names`` to check the order and ``params.index('param_name')`` to get the index of a particular parameter.


Define any derived quantities
-----------------------------

COSMIC requires ``mass_2``, ``separation``, and ``metallicity`` in addition to the
directly-sampled parameters.  The ``compute_derived`` callback converts a ``(N, D)``
array of physical-space samples and a sorted list of parameter names into a dictionary of
``(N,)`` arrays:

.. code-block:: python

    def compute_derived(samples_physical, param_names):
        """Convert sampled parameters to COSMIC inputs."""
        idx = {name: i for i, name in enumerate(param_names)}

        mass_1 = samples_physical[:, idx['mass_1']]
        q      = samples_physical[:, idx['q']]
        porb   = samples_physical[:, idx['porb']]   # physical days after 10^x transform
        z      = samples_physical[:, idx['metallicity']]

        mass_2 = mass_1 * q

        # Kepler's third law: a [AU] from P [yr] and M [M_sun]; then convert to AU
        separation = ((porb / 365.25) ** 2 * (mass_1 + mass_2)) ** (1.0 / 3.0)

        return {
            'mass_2':        mass_2,
            'metallicity_1': z,
            'metallicity_2': z,
            'separation':    separation,   # AU, required by default_reject
        }

.. note::

    The keys ``'mass_2'``, ``'metallicity_1'``, ``'metallicity_2'``, and ``'separation'``
    are expected by :func:`~cosmic.sample.stroopwafel.rejection.default_reject`.  If you
    supply a custom rejection function you are free to use different key names.


Choose a rejection function
---------------------------

Before a batch is passed to ``COSMIC``, unphysical systems are filtered out: stars already
overflowing their Roche lobes at ZAMS, binaries whose components are in contact, and
systems below the hydrogen-burning limit for the secondary.  The built-in
:func:`~cosmic.sample.stroopwafel.rejection.default_reject` function performs all of these
checks:

.. code-block:: python

    from cosmic.sample.stroopwafel.rejection import default_reject

It expects the arrays produced by ``compute_derived`` and returns a boolean mask where
``True`` means the system is rejected.  For most use cases this default is appropriate.
If you need extra cuts — for example, discarding systems with very low metallicity or
imposing a minimum primary mass — you can wrap it:

.. code-block:: python

    def my_reject(samples_physical, derived, param_names):
        base_mask = default_reject(samples_physical, derived, param_names)
        # additionally reject secondaries below 1 M_sun
        base_mask |= (derived['mass_2'] < 1.0)
        return base_mask


Identify what constitutes a hit
-------------------------------

The ``is_interesting`` argument identifies which evolved systems count as hits.  It receives
the ``COSMIC`` ``bpp`` DataFrame for the current batch and must return a tuple
``(n_hits, hit_bin_nums)`` where ``hit_bin_nums`` is an integer array of ``bin_num`` values
(0-indexed within the batch).

The ``STROOPWAFEL`` sampler comes with two preset functions in
:mod:`cosmic.sample.stroopwafel.presets`.

``any_dco(kstar_1, kstar_2)``
    Selects all bound double compact objects whose stellar types match the supplied lists,
    regardless of merger time.  Use this for populations where you care about the DCO
    existing rather than merging within the Hubble time.

``merging_dco(kstar_1, kstar_2, max_merge_time=13.7)``
    Like ``any_dco`` but additionally requires the merger time (computed via LEGWORK) to be
    less than ``max_merge_time`` Gyr.  Suitable for gravitational-wave source studies.

.. code-block:: python

    from cosmic.sample.stroopwafel.presets import any_dco, merging_dco

    # All bound BH-BH systems, no merger time cut
    is_interesting = any_dco(kstar_1=[14], kstar_2=[14])

    # Only BH-BH systems merging within the Hubble time
    is_interesting_merging = merging_dco(kstar_1=[14], kstar_2=[14], max_merge_time=13.7)

See :ref:`kstar-table` for the full list of stellar type codes.

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

With all the pieces in place, you can run the sampler with :class:`~cosmic.sample.stroopwafel.AdaptiveSampler`.  The most important arguments are the parameter space, the total number of systems to evolve, the batch size, the BSE physics settings, the derived quantity function, the rejection function, and the hit function. See the API documentation (:class:`~cosmic.sample.stroopwafel.AdaptiveSampler`) for a full list of options.

The examples below demonstrate how you could go about this.

Examples
--------

Bound BH + BH binaries
^^^^^^^^^^^^^^^^^^^^^^

The following end-to-end example samples all bound BH-BH systems (no merger time
restriction) using a five-dimensional parameter space covering primary mass, mass ratio,
orbital period, eccentricity, and metallicity.

.. include:: ../../../_generated/default_bsedict.rst

.. code-block:: python

    import numpy as np
    from cosmic.sample.stroopwafel import AdaptiveSampler, ParameterSpace, Parameter
    from cosmic.sample.stroopwafel.presets import any_dco
    from cosmic.sample.stroopwafel.rejection import default_reject

    # ------------------------------------------------------------------
    # Parameter space
    # ------------------------------------------------------------------
    params = ParameterSpace([
        Parameter('mass_1',      5.0,    150.0,  sampler='kroupa',     prior='kroupa'),
        Parameter('q',           0.01,   1.0,    sampler='uniform',    prior='uniform'),
        Parameter('porb',        0.15,   5.5,    sampler='sana',       prior='sana'),
        Parameter('ecc',         1e-9,   0.9999, sampler='sana_ecc',   prior='sana_ecc'),
        Parameter('metallicity', 0.0001, 0.03,   sampler='flat_in_log',prior='flat_in_log'),
    ])

    # ------------------------------------------------------------------
    # Derived quantities
    # ------------------------------------------------------------------
    def compute_derived(samples_physical, param_names):
        idx = {name: i for i, name in enumerate(param_names)}
        mass_1 = samples_physical[:, idx['mass_1']]
        q      = samples_physical[:, idx['q']]
        porb   = samples_physical[:, idx['porb']]
        z      = samples_physical[:, idx['metallicity']]
        mass_2     = mass_1 * q
        separation = ((porb / 365.25) ** 2 * (mass_1 + mass_2)) ** (1.0 / 3.0)
        return {
            'mass_2':        mass_2,
            'metallicity_1': z,
            'metallicity_2': z,
            'separation':    separation,
        }

    # ------------------------------------------------------------------
    # Run
    # ------------------------------------------------------------------
    sampler = AdaptiveSampler(
        parameter_space=params,
        total_systems=50_000,
        batch_size=500,
        BSEDict=BSEDict,
        compute_derived=compute_derived,
        reject_systems=default_reject,
        is_interesting=any_dco(kstar_1=[14], kstar_2=[14]),
        output_path='output/bhbh',
        nproc=4,
        n_generations=1,
        seed=42,
    )

    result = sampler.run()

    print(f"Total hits:        {result.num_hits}")
    print(f"Weighted hit rate: {result.hit_rate:.4e} ± {result.hit_rate_uncertainty:.4e}")


BH + star binaries surviving 100 Myr
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For outcomes that are less extreme but still rare — such as persistent BH + star systems —
STROOPWAFEL provides substantial efficiency gains over flat Monte Carlo.  Using the same
parameter space, ``compute_derived``, and ``BSEDict`` as the previous examples, we can simply swap out the hit function to find BH + star systems that remain bound for at least 100 Myr after the BH forms:

.. code-block:: python

    import numpy as np
    from cosmic.sample.stroopwafel import AdaptiveSampler
    from cosmic.sample.stroopwafel.rejection import default_reject

    sampler = AdaptiveSampler(
        parameter_space=params,          # reuse from BHBH example
        total_systems=20_000,
        batch_size=500,
        BSEDict=BSEDict,                 # reuse from BHBH example
        compute_derived=compute_derived, # reuse from BHBH example
        reject_systems=default_reject,
        is_interesting=bh_star_100myr,   # we defined this earlier
        output_path='output/bh_star',
        nproc=4,
        n_generations=1,
        seed=42,
    )

    result = sampler.run()

Because BH + star systems are more common than merging BH-BH pairs, a smaller total budget
is needed and fewer refinement generations are required before the mixture model is
well-constrained.

Rules of thumb
==============

Choosing ``total_systems``, ``batch_size``, and ``n_generations`` is something of an art but there are some rules of thumb to get you started.  The optimal settings depend on the rarity and complexity of the target population, the dimensionality of the parameter space, and your computational resources.

``batch_size``
--------------

``batch_size`` sets how many systems are passed to
:meth:`~cosmic.evolve.Evolve.evolve` per call.

* Aim for ``batch_size`` to be a multiple of ``nproc`` so that COSMIC distributes work
  evenly across cores.
* Values of 200-1000 are typical.  Batches smaller than ~50 increase Python overhead per
  call; batches larger than ~5000 may cause memory pressure on the output DataFrames.

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

Once you have your samples, you can save them to disk as an HDF5 file with the :meth:`~cosmic.output.COSMICSTROOPWAFELResult.save` method.  This saves the parameter samples, derived quantities, and hit information in a compact format that can be loaded later for analysis.

.. code-block:: python

    result.save('bhbh_samples.h5')

We'll talk more about how to load and analyse these results in the :ref:`adaptive_outputs` tutorial next!