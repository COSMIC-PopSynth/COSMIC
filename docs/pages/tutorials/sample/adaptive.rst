.. _adaptive:

**********************************************
Adaptive importance sampling for rare systems
**********************************************

The standard :ref:`independent <independent>` and :ref:`multidimensional <multidim>`
samplers draw binary parameters from their prior distributions and evolve each system with
COSMIC.  For most stellar outcomes this works well: common-envelope episodes, mass-transferring
binaries, and white dwarf systems all occur frequently enough that thousands of random draws
yield a workable sample.

Some outcomes are extremely rare.  Double black holes that merge within the Hubble time
form at rates of order 1-in-10,000 per binary evolved (or far less at high metallicity).
Sampling such populations with a flat Monte Carlo requires tens of millions of COSMIC calls
to accumulate even a few hundred systems — prohibitive for any serious parameter survey.

COSMIC includes a vectorised implementation of the STROOPWAFEL algorithm
(`Broekgaarden et al. 2019 <https://doi.org/10.1093/mnras/stz2558>`_) that solves this
problem using *adaptive importance sampling*.  The sampler first explores parameter space
to locate the progenitor regions of the target population, then concentrates its simulation
budget on those regions, and finally corrects the biased sampling with importance weights so
that any weighted statistic remains an unbiased estimator of the true prior-weighted
distribution.


When should I use this?
========================

Use STROOPWAFEL whenever you need a statistically representative sample of a rare binary
outcome and cannot afford the total binary count that flat Monte Carlo would require.
Typical use cases include:

* Double black holes (or neutron stars) merging within the Hubble time (GW sources)
* BH + stellar-companion systems, such as X-ray binaries or Be/X-ray binaries
* Short-period post-common-envelope binaries that survive to become AM CVn systems
* Any binary channel with a formation efficiency ≲ 10\ :sup:`−3` per prior draw

If your target population is common (≳ 1 % of all binaries evolving as the target), the
plain independent sampler is simpler and fast enough.  The break-even point is roughly
where collecting 100 hits would require more than ~10,000 total evolutions.


How it works: the battleships analogy
======================================

STROOPWAFEL works in three phases that map neatly onto the board game Battleships.

**Exploration — random fire**
    Binaries are drawn at random from the prior distributions and evolved with COSMIC.
    Every binary that produces the desired outcome is recorded as a *hit*.  An adaptive
    stopping criterion (based on the observed hit rate) decides when the remaining budget
    is better spent on refinement than on further random exploration.

**Adaptation — mark the ships**
    One multivariate Gaussian component is placed at each hit location in parameter space.
    The width of each Gaussian is derived from the local density of the prior (via the CDF
    of the prior distribution) so that the proposal is appropriately broad regardless of
    the parameter's scale.  Together the Gaussians form a *mixture model* — a coarse map of
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
    during exploration.  Using these weights, any statistic computed on the hit population
    is an unbiased estimator of the corresponding prior-weighted quantity.


Setting up the parameter space
================================

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

The available samplers and corresponding prior names are:

.. list-table::
    :header-rows: 1
    :widths: 20 45 35

    * - Sampler name
      - Distribution
      - Typical use
    * - ``'uniform'``
      - Uniform between bounds
      - Mass ratio, any flat prior
    * - ``'kroupa'``
      - Kroupa (2001) power law (:math:`dN/dm_1 \propto m_1^{-2.3}`)
      - Primary mass
    * - ``'sana'``
      - Sana et al. (2012) period distribution
        (:math:`dN/d\!\log P \propto (\log P)^{-0.55}`)
      - Orbital period
    * - ``'sana_ecc'``
      - Sana et al. (2012) eccentricity distribution
        (:math:`dN/de \propto e^{-0.45}`)
      - Orbital eccentricity
    * - ``'flat_in_log'``
      - Uniform in :math:`\log_{10}` (Öpik's law)
      - Metallicity, semi-major axis

.. note::

    The ``'sana'`` sampler operates in :math:`\log_{10}(P/\text{days})` space.  The bounds
    you supply are :math:`\log_{10}` values directly — **not** periods in days.  The
    Sana et al. (2012) fit is valid over :math:`\log_{10}(P) \in [0.15,\, 5.5]`,
    corresponding to periods of roughly 1.4 to 316,000 days.

    The lower bound **must be positive** (i.e. :math:`\log_{10}(P_\text{min}) > 0`,
    so :math:`P_\text{min} > 1` day).  Do **not** pass ``np.log10(P_min)`` when that
    value would be negative.

Parameters are stored and returned in alphabetical order by name.  Use
``params.names`` to inspect the column ordering of any sample array.


Defining derived quantities
============================

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


Defining the rejection function
================================

Before a batch is passed to COSMIC, unphysical systems are filtered out: stars already
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


Defining the hit criterion
============================

The ``is_interesting`` argument identifies which evolved systems count as hits.  It receives
the COSMIC ``bpp`` DataFrame for the current batch and must return a tuple
``(n_hits, hit_bin_nums)`` where ``hit_bin_nums`` is an integer array of ``bin_num`` values
(0-indexed within the batch).

STROOPWAFEL ships two preset factory functions in
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

You can also write a fully custom hit function.  For example, to find BH + stellar
companion systems (``kstar_2 ∈ 0–9``) that remain bound for at least 100 Myr after the
BH forms:

.. code-block:: python

    import numpy as np

    _STELLAR_TYPES = set(range(10))   # kstar 0–9: MS through He-giant branch

    def bh_star_100myr(bpp):
        """Hit: BH with a stellar companion bound for at least 100 Myr."""
        bh_star = bpp.loc[
            (
                ((bpp['kstar_1'] == 14) & bpp['kstar_2'].isin(_STELLAR_TYPES))
                | ((bpp['kstar_2'] == 14) & bpp['kstar_1'].isin(_STELLAR_TYPES))
            )
            & (bpp['sep'] > 0)
        ]
        if bh_star.empty:
            return 0, np.array([], dtype=int)

        # Duration in the BH + star state for each binary
        span = bh_star.groupby('bin_num')['tphys'].agg(lambda t: t.max() - t.min())
        hits = span.index[span >= 100.0].values
        return len(hits), hits


Example 1: Bound BH + BH binaries
====================================

The following end-to-end example samples all bound BH-BH systems (no merger time
restriction) using a five-dimensional parameter space covering primary mass, mass ratio,
orbital period, eccentricity, and metallicity.

.. code-block:: python

    import numpy as np
    from cosmic.sample.stroopwafel import AdaptiveSampler, ParameterSpace, Parameter
    from cosmic.sample.stroopwafel.presets import any_dco
    from cosmic.sample.stroopwafel.rejection import default_reject

    # ------------------------------------------------------------------
    # BSE physics settings
    # ------------------------------------------------------------------
    BSEDict = {
        "pts1": 0.001, "pts2": 0.01, "pts3": 0.02, "zsun": 0.014,
        "windflag": 3, "neta": 0.5, "bwind": 0.0, "hewind": 0.5,
        "beta": 0.125, "xi": 0.5, "acc2": 1.5, "LBV_flag": 1,
        "alpha1": 1.0, "lambdaf": 0.0, "ceflag": 1, "cekickflag": 2,
        "cemergeflag": 1, "cehestarflag": 0, "qcflag": 5,
        "qcrit_array": [0.0] * 16,
        "kickflag": 5, "sigma": 265.0, "bhflag": 1, "bhsigmafrac": 1.0,
        "sigmadiv": -20.0, "ecsn": 2.25, "ecsn_mlow": 1.6, "aic": 1,
        "ussn": 1, "polar_kick_angle": 90.0,
        "natal_kick_array": [[-100.0]*5, [-100.0]*5],
        "remnantflag": 4, "mxns": 3.0, "rembar_massloss": 0.5,
        "wd_mass_lim": 1, "grflag": 1, "eddfac": 10, "tflag": 1,
        "ST_tide": 1, "ifflag": 1, "wdflag": 1, "epsnov": 0.001,
        "bdecayfac": 1, "bconst": 3000, "ck": 1000, "htpmb": 1,
        "ST_cr": 1, "rtmsflag": 0,
        "fprimc_array": [2.0 / 21.0] * 16,
        "mm_mu_ns": 400.0, "mm_mu_bh": 200.0, "pisn": -2,
    }

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
        n_generations=3,
        seed=42,
    )

    result = sampler.run()

    print(f"Total hits:        {result.num_hits}")
    print(f"Weighted hit rate: {result.hit_rate:.4e} ± {result.hit_rate_uncertainty:.4e}")


Example 2: BH + star binaries surviving 100 Myr
==================================================

For outcomes that are less extreme but still rare — such as persistent BH + star systems —
STROOPWAFEL provides substantial efficiency gains over flat Monte Carlo.  Using the same
parameter space, ``compute_derived``, and ``BSEDict`` as Example 1:

.. code-block:: python

    import numpy as np
    from cosmic.sample.stroopwafel import AdaptiveSampler
    from cosmic.sample.stroopwafel.rejection import default_reject

    _STELLAR_TYPES = set(range(10))   # kstar 0–9: MS through He-giant

    def bh_star_100myr(bpp):
        """Hit: BH + stellar companion bound for at least 100 Myr."""
        bh_star = bpp.loc[
            (
                ((bpp['kstar_1'] == 14) & bpp['kstar_2'].isin(_STELLAR_TYPES))
                | ((bpp['kstar_2'] == 14) & bpp['kstar_1'].isin(_STELLAR_TYPES))
            )
            & (bpp['sep'] > 0)
        ]
        if bh_star.empty:
            return 0, np.array([], dtype=int)
        span = bh_star.groupby('bin_num')['tphys'].agg(lambda t: t.max() - t.min())
        hits = span.index[span >= 100.0].values
        return len(hits), hits

    sampler = AdaptiveSampler(
        parameter_space=params,          # reuse from Example 1
        total_systems=20_000,
        batch_size=500,
        BSEDict=BSEDict,                 # reuse from Example 1
        compute_derived=compute_derived, # reuse from Example 1
        reject_systems=default_reject,
        is_interesting=bh_star_100myr,
        output_path='output/bh_star',
        nproc=4,
        n_generations=2,
        seed=42,
    )

    result = sampler.run()

Because BH + star systems are more common than merging BH-BH pairs, a smaller total budget
is needed and fewer refinement generations are required before the mixture model is
well-constrained.


Choosing ``total_systems``, ``batch_size``, and ``n_generations``
===================================================================

``batch_size``
--------------

``batch_size`` sets how many systems are passed to
:meth:`~cosmic.evolve.Evolve.evolve` per call.

* Aim for ``batch_size`` to be a multiple of ``nproc`` so that COSMIC distributes work
  evenly across cores.
* Values of 200–1000 are typical.  Batches smaller than ~50 increase Python overhead per
  call; batches larger than ~5000 may cause memory pressure on the output DataFrames.
* A practical starting point is ``batch_size = 100 * nproc``.

``total_systems``
-----------------

This is the total number of binary evolutions across all phases.

* For **very rare events** (hit rate ≲ 10\ :sup:`−4`, e.g. merging BH-BH at near-solar
  metallicity), start with ``total_systems`` in the range 100,000–500,000.  The exploration
  phase will find tens to hundreds of hits; refinement then multiplies that count many-fold.
* For **moderately rare events** (hit rate ~ 10\ :sup:`−3` to 10\ :sup:`−2`, e.g. any bound
  BH-BH or long-lived BH + star), 20,000–50,000 systems is usually sufficient.
* As a rule of thumb, aim for at least ~30 hits during exploration before the adaptation
  phase begins — fewer hits lead to a poorly-constrained Gaussian mixture.  If exploration
  ends with very few hits, increase ``total_systems`` and re-run.

``n_generations``
-----------------

Each refinement generation uses an equal share of the remaining budget after exploration.
The EM step between generations can improve the mixture, but with diminishing returns:

* ``n_generations = 1`` (the default) uses the mixture as constructed from exploration
  hits, with no EM updates.  This is a good starting point for any new target population.
* ``n_generations = 3`` gives a noticeable improvement for very rare populations with
  complex progenitor structure.
* Beyond 5 generations the returns diminish rapidly.

.. tip::

    To run a plain Monte Carlo without any adaptation or refinement — useful as a baseline
    or for common populations — pass ``mc_only=True``.  The sampler will draw from the
    prior only, and importance weights will equal the prior density divided by itself
    (i.e. all weights are equal).


Working with results
======================

:meth:`~cosmic.sample.stroopwafel.engine.AdaptiveSampler.run` returns a
:class:`~cosmic.sample.stroopwafel.result.STROOPWAFELResult` object containing all
simulated systems, their importance weights, and summary statistics:

.. code-block:: python

    import numpy as np

    # Shape of sample array and column ordering
    print(result.samples.shape)   # (N_total, D)
    print(result.param_names)     # sorted alphabetically, e.g.
                                  # ['ecc', 'mass_1', 'metallicity', 'porb', 'q']

    # Extract hits and their weights
    hit_samples = result.samples[result.is_hit]        # (N_hits, D)
    hit_weights = result.weights[result.is_hit]        # (N_hits,)

    # Normalise weights for the hit population
    hit_weights_norm = hit_weights / hit_weights.sum()

    # Importance-weighted primary mass histogram
    m1_col = result.param_names.index('mass_1')
    m1_hits = hit_samples[:, m1_col]
    hist, edges = np.histogram(m1_hits, bins=20, weights=hit_weights_norm)

    # Importance-weighted hit rate (fraction of prior draws producing a hit)
    print(f"Hit rate: {result.hit_rate:.4e} ± {result.hit_rate_uncertainty:.4e}")

    # How the budget was spent
    print(f"Explored: {result.num_explored}   Total hits: {result.num_hits}")

.. note::

    ``result.samples`` stores samples in **physical space** (masses in M\ :sub:`☉`,
    periods in days, etc.) in the alphabetically-sorted column order defined by
    ``ParameterSpace``.  Always use ``result.param_names`` to map column indices to
    parameter names rather than relying on the order in which you defined the parameters.


Saving and loading results
============================

The full result can be saved to HDF5 for later analysis:

.. code-block:: python

    from cosmic.sample.stroopwafel import io as swio

    swio.save_result('bhbh_result.h5', result)

The file stores the sample array, importance weights, hit flags, generation labels, and
summary statistics as HDF5 datasets and attributes.

To reload the data in a later session without re-running the sampler:

.. code-block:: python

    import h5py
    import numpy as np

    with h5py.File('bhbh_result.h5', 'r') as f:
        samples     = f['samples'][:]
        weights     = f['weights'][:]
        is_hit      = f['is_hit'][:]
        param_names = list(f.attrs['param_names'])
        num_hits    = f.attrs['num_hits']

    hits = samples[is_hit]
    m1   = hits[:, param_names.index('mass_1')]
    w    = weights[is_hit]
    print(f"Weighted mean BH primary mass: {np.average(m1, weights=w):.1f} M_sun")

.. note::

    The HDF5 file does **not** store the raw COSMIC ``bpp`` output tables.  If you need
    the evolutionary histories of the hit systems, re-evolve them with COSMIC using the
    hit sample coordinates from ``result.samples[result.is_hit]`` and the same ``BSEDict``.
