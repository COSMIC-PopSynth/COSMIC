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
after adaptation — once the Gaussian mixture has been fit to the exploration hits — and
resume the (usually much larger) refinement phase separately. This is useful when you want
to

* run exploration and refinement as **separate cluster jobs**, perhaps with different
  walltimes or allocations;
* spread a single adapted mixture out across **several refinement jobs** on different nodes; or
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

    params = ParameterSpace([
        Parameter('mass_1',      5.0,        150.0,      dist='kroupa'),
        Parameter('q',           0.01,       1.0,        dist='uniform'),
        Parameter('porb',        10**(0.15), 10**(5.5),  dist='sana'),
        Parameter('ecc',         1e-9,       0.9999,     dist='sana_ecc'),
        Parameter('metallicity', 0.0001,     0.03,       dist='flat_in_log'),
    ])

    def derive_params(sampled):
        return {'mass_2': sampled['mass_1'] * sampled['q']}

    sampler = AdaptiveSampler(
        parameter_space=params,
        total_systems=500_000,
        batch_size=1000,
        BSEDict=BSEDict,
        SSEDict={'stellar_engine': 'sse'},
        is_interesting=any_dco(kstar_1=[14], kstar_2=[14]),
        derive_params=derive_params,
        reject_systems="default",
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

    sampler = AdaptiveSampler.from_checkpoint('checkpoint.h5')
    result = sampler.run_refinement()
    result.save('result.h5')

You do not need to re-import or re-specify the parameter space, ``BSEDict``, ``SSEDict``, or
any of the callables — they were all saved into the checkpoint.  If you *want* to change something for
the refinement phase (a common one is running on more cores, or with a larger budget than
exploration), pass it as a keyword override:

.. code-block:: python

    sampler = AdaptiveSampler.from_checkpoint(
        'checkpoint.h5', nproc=16
    )

Any :class:`~cosmic.sample.stroopwafel.AdaptiveSampler` constructor argument may be
overridden this way (``parameter_space``, ``BSEDict``, ``SSEDict``, ``derive_params``,
``reject_systems``, ``is_interesting``, ``batch_size``, ``nproc``,
``kappa``, ``n_generations``, ``only_save_hit_tables``, ``seed``).

The ``result`` is an ordinary :class:`~cosmic.output.COSMICStroopOutput` — identical in form
to what a single :meth:`~cosmic.sample.stroopwafel.AdaptiveSampler.run` would have produced —
so you can analyse it exactly as described in :ref:`adaptive_outputs`.


What is stored in a checkpoint
==============================

A checkpoint is a complete snapshot — it stores both the exploration *results* and the full
*configuration*, so refinement can resume with no further input:

* the fitted Gaussian mixture model;
* every exploration sample (in the internal sampling space) and its hit/bookkeeping flags;
* the COSMIC output tables (``bpp``, ``bcm``, ``initC``, ``kick_info``) for the explored
  systems;
* the scalar counters needed to compute unbiased weights later (``num_explored``,
  ``fraction_explored``, ``prior_fraction_rejected``, ...); and
* the full configuration needed to rebuild the sampler: the parameter space, ``BSEDict``,
  ``SSEDict``, the ``derive_params`` / ``reject_systems`` / ``is_interesting`` callables, the
  remaining scalar settings, and the live RNG state.

The callables and parameter space are serialised with :mod:`dill` so it lets you use general functions.
Because the RNG state is stored too, refinement continues the random stream seamlessly rather than restarting it (pass ``seed=`` to
``from_checkpoint`` if you instead want a fresh stream).


Reusing one checkpoint for several refinement jobs
==================================================

Because :meth:`~cosmic.sample.stroopwafel.AdaptiveSampler.from_checkpoint` reads the
checkpoint without modifying it, you can launch any number of independent refinement jobs
from the same ``checkpoint.h5`` and each will draw a fresh, independent refinement sample from the shared mixture.
You can also override the budget stored in the checkpoint by passing ``total_systems`` or
``n_generations`` to :meth:`~cosmic.sample.stroopwafel.AdaptiveSampler.from_checkpoint`.

.. warning::

    You should make sure that you change the ``seed`` when you launch multiple refinement jobs from the same checkpoint, otherwise they will all draw the same random numbers and produce identical results. The simplest way to do this is to pass ``seed=None`` to :meth:`~cosmic.sample.stroopwafel.AdaptiveSampler.from_checkpoint`, which will seed the RNG from the system clock.


Wrap-up
=======

And that's all on adaptive sampling folks! You should now know everything you need to run an adaptive importance sampling simulation in ``COSMIC`` - enjoy exploring those rare populations!