.. _adaptive_distributions:

*********************************
Distributions and custom priors
*********************************

This tutorial assumes that you've already gone through :ref:`adaptive_basics`.

When using adaptive sampling you need to define a distribution to use for each :class:`~cosmic.sample.stroopwafel.Parameter` that you sample.
The distribution performs three operations: it draws samples, it evaluates the prior probability density, and it defines
the adaptive-sampling kernel width used during refinement.

In this tutorial we'll cover the built-in distributions and show how to define your own custom distributions and transforms.

Built-in distributions
=======================

You can select any of the following distributions as the ``dist`` argument to a :class:`~cosmic.sample.stroopwafel.Parameter`:

.. list-table::
    :header-rows: 1
    :widths: 18 22 60

    * - Name
      - Sampled as
      - Description
    * - ``'uniform'``
      - uniform
      - Flat between the bounds.  Good default for mass ratio, etc.
    * - ``'flat_in_log'``
      - uniform in :math:`\log_{10}`
      - Flat in the log of the parameter (e.g. metallicity).
    * - ``'kroupa'``
      - broken power law (:math:`\alpha = -1.3` for :math:`m < 0.5\,M_\odot`,
        :math:`-2.3` above)
      - Kroupa initial mass function for the primary mass.  Use a lower bound of
        :math:`\geq 0.08\,M_\odot` (COSMIC cannot evolve lower-mass stars).
    * - ``'sana'``
      - power law in :math:`\log_{10} P`, :math:`\alpha = -0.55`
      - Sana et al. (2012) orbital-period distribution.
    * - ``'sana_ecc'``
      - power law, :math:`\alpha = -0.45`
      - Sana et al. (2012) eccentricity distribution.
    * - ``'uniform_in_sine'``
      - uniform in :math:`\sin\theta`
      - Isotropic angle (e.g. inclination-like coordinates).
    * - ``'uniform_in_cosine'``
      - uniform in :math:`\cos\theta`
      - Isotropic angle for declination-like coordinates.
    * - ``'disberg'``
      - log-normal
      - Natal-kick magnitude, :math:`\ln v \sim \mathcal{N}(5.67, 0.59)`.

If you ever want to access this, you can get the full list of registered distributions with :data:`cosmic.sample.stroopwafel.distributions.DISTRIBUTIONS`.


How distributions are built
===========================

Under the hood, we set each distribution up as a combination of a base distribution and a coordinate transform. This allows you to mix and match options:

* the base distribution (:class:`~cosmic.sample.stroopwafel.distributions.Uniform`,
  :class:`~cosmic.sample.stroopwafel.distributions.PowerLaw`,
  :class:`~cosmic.sample.stroopwafel.distributions.BrokenPowerLaw`, or
  :class:`~cosmic.sample.stroopwafel.distributions.TruncatedNormal`) handles sampling and
  the density in the *sampling space*; and
* the transform (:class:`~cosmic.sample.stroopwafel.distributions.Log10`,
  :class:`~cosmic.sample.stroopwafel.distributions.Ln`, etc.) maps between the physical
  space you specify bounds in and that sampling space.

This is why, for example, ``'flat_in_log'`` is just a uniform distribution paired with a
:math:`\log_{10}` transform, and ``'sana'`` is a power law paired with the same transform.
You can build the same objects yourself:

.. code-block:: python

    from cosmic.sample.stroopwafel.distributions import (
        Uniform, PowerLaw, BrokenPowerLaw, TruncatedNormal, Log10,
    )

    Uniform(transform=Log10())                         # equivalent to 'flat_in_log'
    PowerLaw(-0.55, transform=Log10())                 # equivalent to 'sana'
    BrokenPowerLaw(breaks=[0.5], alphas=[-1.3, -2.3])  # equivalent to 'kroupa'

Bounds are always given to a :class:`~cosmic.sample.stroopwafel.Parameter` in **physical**
space; the transform converts them into sampling space automatically.


Defining your own distribution
==============================

Now let's say that you want to define your own distribution. You can do this in three ways - let's take a look at them in order of increasing complexity.

1. Pass a distribution instance directly
----------------------------------------

The quickest option is to tweak one of the built-in base distributions like we did above and hand the
instance straight to a :class:`~cosmic.sample.stroopwafel.Parameter` via ``dist``.

.. code-block:: python

    from cosmic.sample.stroopwafel import Parameter
    from cosmic.sample.stroopwafel.distributions import PowerLaw

    # a steeper-than-Kroupa IMF for the primary mass
    Parameter('mass_1', 5.0, 150.0, dist=PowerLaw(-2.7))


2. Register a named distribution
--------------------------------

If you want to reuse a distribution across several parameter spaces — or simply refer to it
by a memorable name — you can register it once with
:func:`~cosmic.sample.stroopwafel.distributions.register`. This then allows you to refer to it by name in any :class:`~cosmic.sample.stroopwafel.Parameter`:

.. code-block:: python

    from cosmic.sample.stroopwafel import Parameter
    from cosmic.sample.stroopwafel.distributions import register, PowerLaw

    register('imf_steep', PowerLaw(-2.7))

    # ... anywhere later
    Parameter('mass_1', 5.0, 150.0, dist='imf_steep')


3. Write a new distribution class
---------------------------------

For a genuinely new functional form, you'll need to create a new class that subclasses off
:class:`~cosmic.sample.stroopwafel.distributions.Distribution` and implement two methods:

``sample(n, lo, hi, rng)``
    Draw ``n`` samples in **sampling space**, restricted to ``[lo, hi]``.

``pdf(values, lo, hi)``
    Return the prior density at ``values``, normalised over ``[lo, hi]`` in sampling space.

Both ``lo`` and ``hi`` are bounds in sampling space — the parameter space has already
applied the transform, so you do not need to worry about it here.  The example below
implements a truncated exponential distribution:

.. code-block:: python

    import numpy as np
    from cosmic.sample.stroopwafel import Parameter
    from cosmic.sample.stroopwafel.distributions import Distribution

    class Exponential(Distribution):
        """p(x) ∝ exp(-x / scale), truncated to [lo, hi]."""

        def __init__(self, scale, transform=None):
            super().__init__(transform)
            self.scale = scale

        def sample(self, n, lo, hi, rng=None):
            rng = rng or np.random.default_rng()
            u = rng.uniform(0, 1, n)
            c_lo, c_hi = np.exp(-lo / self.scale), np.exp(-hi / self.scale)
            return -self.scale * np.log(c_lo - u * (c_lo - c_hi))   # inverse CDF

        def pdf(self, values, lo, hi):
            norm = self.scale * (np.exp(-lo / self.scale) - np.exp(-hi / self.scale))
            return np.exp(-values / self.scale) / norm

    Parameter('some_param', 0.0, 10.0, dist=Exponential(scale=2.0))

And this class defines everything we need, our distribution now works everywhere the built-ins do, and
can be combined with any transform (``dist=Exponential(2.0, transform=Log10())``).

.. note::

    During refinement, ``STROOPWAFEL`` places a Gaussian kernel at each hit whose width is
    set by :meth:`~cosmic.sample.stroopwafel.distributions.Distribution.sigma`.  The default
    implementation, ``avg_density / pdf``, is appropriate for almost all distributions and
    is inherited automatically — you only need to override it if you have an exact
    closed-form CDF you would rather step through (as
    :class:`~cosmic.sample.stroopwafel.distributions.PowerLaw` does).


Custom transforms
=================

Transforms are just as extensible.  A transform implements ``to_sampling`` (physical →
sampling) and ``to_physical`` (sampling → physical); the corresponding bound conversion is
derived automatically and handles decreasing maps by swapping the endpoints.  The built-in
transforms are :class:`~cosmic.sample.stroopwafel.distributions.Identity`,
:class:`~cosmic.sample.stroopwafel.distributions.Log10`,
:class:`~cosmic.sample.stroopwafel.distributions.Ln`,
:class:`~cosmic.sample.stroopwafel.distributions.Sin`, and
:class:`~cosmic.sample.stroopwafel.distributions.CosShift`.  To add your own, subclass
:class:`~cosmic.sample.stroopwafel.distributions.Transform`:

.. code-block:: python

    import numpy as np
    from cosmic.sample.stroopwafel.distributions import Transform

    class Sqrt(Transform):
        def to_sampling(self, values):
            return np.sqrt(values)

        def to_physical(self, values):
            return values ** 2

Wrap-up
=======

And that's everything you need to know about distributions and transforms in ``COSMIC``'s implementation of '``STROOPWAFEL``.
You can now define your own custom priors and use them in your adaptive sampling runs.

Next, we'll look at how to analyse the outputs of an adaptive sampling run in :ref:`adaptive_outputs`.