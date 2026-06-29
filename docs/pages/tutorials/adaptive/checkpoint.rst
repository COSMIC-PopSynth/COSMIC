.. _adaptive_checkpoint:

******************************
Saving and loading checkpoints
******************************

In this tutorial, we will cover how to save and load checkpoints during an adaptive importance sampling run in ``COSMIC``. This allows you to save the state of your simulation once the adaptation phase is complete, and then load it later to continue sampling from the same Gaussian mixture model (GMM) without having to redo the adaptation phase. This can be particularly useful if you want to run multiple sampling phases with the same adapted GMM (e.g. over multiple nodes on a cluster).

This tutorial assumes that you've already gone through :ref:`adaptive_basics`.

Why split a run in two?
=======================

A call to :meth:`~cosmic.sample.stroopwafel.AdaptiveSampler.run` performs all three phases
back to back: exploration, adaptation, and refinement.  Splitting the run lets you stop
after adaptation — once the Gaussian mixture has been fitted to the exploration hits — and
resume the (usually much larger) refinement phase separately.  This is useful when you want
to

* run exploration and refinement as **separate cluster jobs**, perhaps with different
  walltimes or allocations;
* fan a single adapted mixture out across **several refinement jobs** on different nodes; or
* simply **inspect** the mixture before committing compute to refinement.

Instead of :meth:`~cosmic.sample.stroopwafel.AdaptiveSampler.run`, you use the two
multi-job entry points:
:meth:`~cosmic.sample.stroopwafel.AdaptiveSampler.run_exploration` (which returns a
checkpoint) and :meth:`~cosmic.sample.stroopwafel.AdaptiveSampler.run_refinement`.


Stage 1 — explore, adapt, and save a checkpoint
===============================================

Set up the sampler exactly as you would for a normal run (see :ref:`adaptive_basics`), but
call :meth:`~cosmic.sample.stroopwafel.AdaptiveSampler.run_exploration` instead of
``run()``.  This runs the exploration and adaptation phases and returns a
:class:`~cosmic.output.STROOPWAFELCheckpoint`, which you then save to disk.

.. code-block:: python

    from cosmic.sample.stroopwafel import AdaptiveSampler, ParameterSpace, Parameter
    from cosmic.sample.stroopwafel.presets import any_dco
    from cosmic.sample.stroopwafel.rejection import default_reject

    params = ParameterSpace([
        Parameter('mass_1',      5.0,        150.0,      dist='kroupa'),
        Parameter('q',           0.01,       1.0,        dist='uniform'),
        Parameter('porb',        10**(0.15), 10**(5.5),  dist='sana'),
        Parameter('ecc',         1e-9,       0.9999,     dist='sana_ecc'),
        Parameter('metallicity', 0.0001,     0.03,       dist='flat_in_log'),
    ])

    def compute_derived(samples_physical, param_names):
        idx = {name: i for i, name in enumerate(param_names)}
        mass_1 = samples_physical[:, idx['mass_1']]
        q      = samples_physical[:, idx['q']]
        porb   = samples_physical[:, idx['porb']]
        z      = samples_physical[:, idx['metallicity']]
        mass_2     = mass_1 * q
        separation = ((porb / 365.25) ** 2 * (mass_1 + mass_2)) ** (1.0 / 3.0)
        return {'mass_2': mass_2, 'metallicity_1': z,
                'metallicity_2': z, 'separation': separation}

    sampler = AdaptiveSampler(
        parameter_space=params,
        total_systems=500_000,
        batch_size=1000,
        BSEDict=BSEDict,
        compute_derived=compute_derived,
        reject_systems=default_reject,
        is_interesting=any_dco(kstar_1=[14], kstar_2=[14]),
        output_path='output/explore',
        nproc=4,
        seed=42,
    )

    checkpoint = sampler.run_exploration()   # exploration + adaptation only
    checkpoint.save('checkpoint.h5')

.. note::

    If exploration finds no hits there is nothing to adapt, and the checkpoint's mixture
    will be ``None``.  In that case you should increase ``total_systems`` and re-run
    exploration before attempting refinement (see the rules of thumb in
    :ref:`adaptive_basics`).


Stage 2 — load the checkpoint and refine
========================================

In a second script (or cluster job) rebuild the sampler with
:meth:`~cosmic.sample.stroopwafel.AdaptiveSampler.from_checkpoint`, then call
:meth:`~cosmic.sample.stroopwafel.AdaptiveSampler.run_refinement`.

.. code-block:: python

    from cosmic.sample.stroopwafel import AdaptiveSampler
    from cosmic.sample.stroopwafel.presets import any_dco
    from cosmic.sample.stroopwafel.rejection import default_reject

    # `params`, `compute_derived`, and `BSEDict` must be available again here —
    # in practice, import them from a shared module used by both jobs.

    sampler = AdaptiveSampler.from_checkpoint(
        'checkpoint.h5',
        parameter_space=params,
        batch_size=1000,
        BSEDict=BSEDict,
        compute_derived=compute_derived,
        reject_systems=default_reject,
        is_interesting=any_dco(kstar_1=[14], kstar_2=[14]),
        output_path='output/refine',
        nproc=4,
        seed=7,
    )

    result = sampler.run_refinement()
    result.save('result.h5')

The ``result`` is an ordinary :class:`~cosmic.output.COSMICStroopOutput` — identical in form
to what a single :meth:`~cosmic.sample.stroopwafel.AdaptiveSampler.run` would have produced —
so you can analyse it exactly as described in :ref:`adaptive_outputs`.


What is and isn't stored in a checkpoint
========================================

A checkpoint stores everything that is *derived from running COSMIC*:

* the fitted Gaussian mixture model,
* every exploration sample (in the internal sampling space) and its hit/bookkeeping flags,
* the COSMIC output tables (``bpp``, ``bcm``, ``initC``, ``kick_info``) for the explored
  systems, and
* the scalar counters needed to compute unbiased weights later (``num_explored``,
  ``fraction_explored``, ``prior_fraction_rejected``, ``total_systems``, ...).

It deliberately does **not** store your Python callables or physics settings — the
``BSEDict``, ``compute_derived``, ``reject_systems``, and ``is_interesting`` functions are
not serialisable in general, so you must supply them again to
:meth:`~cosmic.sample.stroopwafel.AdaptiveSampler.from_checkpoint`.  Only the parameter
**names** are stored, and they are checked against the ``parameter_space`` you provide; a
mismatch raises a ``ValueError`` to stop you from refining against the wrong setup.


Reusing one checkpoint for several refinement jobs
==================================================

Because :meth:`~cosmic.sample.stroopwafel.AdaptiveSampler.from_checkpoint` reads the
checkpoint without modifying it, you can launch any number of independent refinement jobs
from the same ``checkpoint.h5`` — for example with different ``seed`` values on different
nodes — and each will draw a fresh, independent refinement sample from the shared mixture.
You can also override the budget stored in the checkpoint by passing ``total_systems`` or
``n_generations`` to :meth:`~cosmic.sample.stroopwafel.AdaptiveSampler.from_checkpoint`.

