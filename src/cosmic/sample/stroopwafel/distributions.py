"""Composable distributions for STROOPWAFEL parameter sampling.

A :class:`Distribution` bundles the four things STROOPWAFEL needs to know
about a parameter's prior, which used to be spread across ``samplers.py``,
``priors.py`` and ``transforms.py``:

* how to **draw** samples (in sampling space),
* the prior **pdf** (in sampling space),
* the adaptive-sampling **sigma** (Gaussian kernel width), and
* the **transform** between physical and sampling space.

Each distribution is a *base distribution* (:class:`Uniform`,
:class:`PowerLaw`, :class:`BrokenPowerLaw`, :class:`TruncatedNormal`) composed
with a *transform* (:class:`Identity`, :class:`Log10`, :class:`Ln`,
:class:`Sin`, :class:`CosShift`).  The base operates entirely in sampling
space; the
transform maps to and from the physical space the user specifies bounds in
and the simulation consumes.  For example ``flat_in_log`` is
``Uniform(transform=Log10())`` and ``sana`` is
``PowerLaw(SANA_G, transform=Log10())``.

To add a distribution, build an instance and either pass it straight to a
:class:`~cosmic.sample.stroopwafel.parameter_space.Parameter` or register it
under a name::

    from cosmic.sample.stroopwafel.distributions import PowerLaw, register
    register('my_imf', PowerLaw(-2.7))

All array-valued methods take and return (N,) ndarrays.
"""
import numpy as np
from scipy.stats import norm as _scipy_norm


# ---------------------------------------------------------------------------
# Transforms between physical and sampling space
# ---------------------------------------------------------------------------
class Transform:
    """Map between physical space and sampling space.

    The base class is the identity transform; subclasses override
    :meth:`to_sampling` and :meth:`to_physical`.  ``bounds`` is derived
    automatically and handles non-monotonic-increasing maps (e.g.
    :class:`CosShift`) by sorting the endpoints.
    """

    def to_sampling(self, values):
        """Convert physical-space values to sampling space.

        Parameters
        ----------
        values : `numpy.ndarray` or `float`
            Value(s) in physical space.

        Returns
        -------
        `numpy.ndarray` or `float`
            Value(s) in sampling space.
        """
        return values

    def to_physical(self, values):
        """Convert sampling-space values to physical space.

        Parameters
        ----------
        values : `numpy.ndarray` or `float`
            Value(s) in sampling space.

        Returns
        -------
        `numpy.ndarray` or `float`
            Value(s) in physical space.
        """
        return values

    def bounds(self, lo, hi):
        """Map a physical-space bound pair into sampling space.

        Parameters
        ----------
        lo : `float`
            Lower bound in physical space.
        hi : `float`
            Upper bound in physical space.

        Returns
        -------
        lo_sampling : `float`
            Lower bound in sampling space.
        hi_sampling : `float`
            Upper bound in sampling space.
        """
        a = self.to_sampling(lo)
        b = self.to_sampling(hi)
        return (a, b) if a <= b else (b, a)


class Identity(Transform):
    """Identity transform: sampling space is physical space."""


class Log10(Transform):
    """Sampling space is ``log10(physical)`` (physical must be positive)."""

    def to_sampling(self, values):
        return np.log10(values)

    def to_physical(self, values):
        return np.power(10.0, values)


class Ln(Transform):
    """Sampling space is ``ln(physical)`` (physical must be positive)."""

    def to_sampling(self, values):
        return np.log(values)

    def to_physical(self, values):
        return np.exp(values)


class Sin(Transform):
    """Sampling space is ``sin(angle)`` for an angle in ``[-pi/2, pi/2]``."""

    def to_sampling(self, values):
        return np.sin(values)

    def to_physical(self, values):
        return np.arcsin(values)


class CosShift(Transform):
    """Sampling space is ``cos(angle + pi/2) = -sin(angle)``.

    Suitable for declination-like coordinates measured from ``-pi/2`` to
    ``pi/2``.  The map is monotonically decreasing, so :meth:`Transform.bounds`
    swaps the endpoints to keep ``lo <= hi``.
    """

    def to_sampling(self, values):
        return np.cos(values + np.pi / 2)

    def to_physical(self, values):
        return np.arccos(values) - np.pi / 2


# ---------------------------------------------------------------------------
# Base distribution
# ---------------------------------------------------------------------------
class Distribution:
    """A prior distribution, sampled and evaluated in sampling space.

    Subclasses implement :meth:`sample` and :meth:`pdf`.  :meth:`sigma` has a
    default implementation (``avg_density / pdf``) that power laws override.
    All ``lo``/``hi`` arguments are bounds in *sampling* space (i.e. already
    passed through ``self.transform``); :class:`ParameterSpace` handles that
    conversion via :meth:`Transform.bounds`.

    Parameters
    ----------
    transform : `Transform`, optional
        Map between physical and sampling space, by default :class:`Identity`.
    """

    def __init__(self, transform=None):
        self.transform = transform if transform is not None else Identity()

    def sample(self, n, lo, hi, rng=None):
        """Draw ``n`` samples in sampling space within ``[lo, hi]``.

        Parameters
        ----------
        n : `int`
            Number of samples to draw.
        lo : `float`
            Lower bound in sampling space.
        hi : `float`
            Upper bound in sampling space.
        rng : `numpy.random.Generator`, optional
            Random number generator, by default None.

        Returns
        -------
        `numpy.ndarray`
            (N,) array of samples in sampling space.
        """
        raise NotImplementedError

    def pdf(self, values, lo, hi):
        """Prior probability density at ``values`` (sampling space).

        Parameters
        ----------
        values : `numpy.ndarray`
            (N,) array of values in sampling space.
        lo : `float`
            Lower bound in sampling space.
        hi : `float`
            Upper bound in sampling space.

        Returns
        -------
        `numpy.ndarray`
            (N,) array of prior densities, normalised over ``[lo, hi]``.
        """
        raise NotImplementedError

    def sigma(self, values, lo, hi, avg_density):
        """Per-sample Gaussian kernel width for adaptive refinement.

        The default places a kernel whose width is the local inter-sample
        spacing, ``avg_density / pdf``.  Distributions with closed-form CDFs
        (e.g. :class:`PowerLaw`) override this with an exact CDF-space step.

        Parameters
        ----------
        values : `numpy.ndarray`
            (K,) array of hit values in sampling space.
        lo : `float`
            Lower bound in sampling space.
        hi : `float`
            Upper bound in sampling space.
        avg_density : `float`
            Characteristic inter-sample spacing.

        Returns
        -------
        `numpy.ndarray`
            (K,) array of sigma values.
        """
        return avg_density / self.pdf(values, lo, hi)


# ---------------------------------------------------------------------------
# Concrete distributions
# ---------------------------------------------------------------------------
class Uniform(Distribution):
    """Uniform distribution on ``[lo, hi]`` in sampling space.

    Combined with a transform this covers ``uniform`` (identity),
    ``flat_in_log`` (:class:`Log10`), ``uniform_in_sine`` (:class:`Sin`) and
    ``uniform_in_cosine`` (:class:`CosShift`).
    """

    def sample(self, n, lo, hi, rng=None):
        rng = rng or np.random.default_rng()
        return rng.uniform(lo, hi, n)

    def pdf(self, values, lo, hi):
        return np.full(len(values), 1.0 / (hi - lo))


class PowerLaw(Distribution):
    r"""Power-law distribution ``p(x) \propto x^alpha`` on ``[lo, hi]``.

    Sampling uses the inverse CDF and the prior is the normalised power law,
    both in sampling space.  Combined with a transform this covers ``kroupa``
    and ``sana_ecc`` (identity) and ``sana`` (:class:`Log10`).

    Parameters
    ----------
    alpha : `float`
        Power-law exponent.
    transform : `Transform`, optional
        Map between physical and sampling space, by default :class:`Identity`.
    """

    def __init__(self, alpha, transform=None):
        super().__init__(transform)
        self.alpha = alpha

    def sample(self, n, lo, hi, rng=None):
        rng = rng or np.random.default_rng()
        u = rng.uniform(0, 1, n)
        a = self.alpha + 1
        return np.power(u * (hi**a - lo**a) + lo**a, 1.0 / a)

    def pdf(self, values, lo, hi):
        a = self.alpha
        norm = (a + 1) / (hi**(a + 1) - lo**(a + 1))
        return norm * np.power(values, a)

    def sigma(self, values, lo, hi, avg_density):
        """Exact CDF-space step for the power-law sigma.

        Maps each hit to its CDF position, steps by ``avg_density`` in CDF
        space, maps back, and returns the larger of the two distances.  The
        CDF of ``p(x) \\propto x^alpha`` on ``[lo, hi]`` is
        ``F(x) = (x^a - lo^a) / (hi^a - lo^a)`` with ``a = alpha + 1``, so the
        inverse is ``F^{-1}(u) = (u*(hi^a - lo^a) + lo^a)^{1/a}``.  Keeping
        every intermediate in ``[lo^a, hi^a]`` avoids the catastrophic
        cancellation seen with very small ``lo``.

        Raises
        ------
        `ValueError`
            If ``lo <= 0`` (the power-law variable must be positive over the
            whole range).
        """
        if lo <= 0:
            raise ValueError(
                f"PowerLaw.sigma requires a positive lower bound in sampling "
                f"space, but got lo={lo}.  When combined with a log transform "
                f"(e.g. 'sana'), make sure the physical lower bound maps to a "
                f"positive sampling-space value (for 'sana', period > 1 day)."
            )
        a = self.alpha + 1
        lo_a = float(lo) ** a
        hi_a = float(hi) ** a
        range_a = hi_a - lo_a  # always finite; no huge intermediates

        # Forward CDF, step in CDF space, then inverse CDF.
        u = np.clip((np.power(values, a) - lo_a) / range_a, 0.0, 1.0)
        u_right = np.clip(u + avg_density, 0.0, 1.0)
        u_left = np.clip(u - avg_density, 0.0, 1.0)
        x_right = np.power(u_right * range_a + lo_a, 1.0 / a)
        x_left = np.power(u_left * range_a + lo_a, 1.0 / a)

        return np.maximum(np.abs(x_right - values), np.abs(x_left - values))


class BrokenPowerLaw(Distribution):
    r"""Continuous broken power law: ``p(x) \propto x^alpha`` with ``alpha``
    changing at fixed breakpoints.

    The density is continuous across every breakpoint.  With the default
    :class:`Identity` transform this gives the Kroupa IMF
    (``breaks=[0.5]``, ``alphas=[-1.3, -2.3]``): shallower below 0.5 Msun,
    steeper above.  Over any window that contains no breakpoint it reduces
    exactly to :class:`PowerLaw`, so sampling masses above 0.5 Msun behaves
    identically to a single power law.

    Sampling and ``sigma`` use a piecewise inverse CDF; the normalisation,
    sampling, and density are all computed over the requested ``[lo, hi]``
    window only (segments outside it are ignored).

    Parameters
    ----------
    breaks : sequence of `float`
        Internal breakpoints in sampling space, strictly increasing.  The
        distribution has ``len(breaks) + 1`` segments.
    alphas : sequence of `float`
        Power-law exponent for each segment (``len(breaks) + 1`` of them);
        ``alphas[i]`` applies below ``breaks[i]`` and ``alphas[-1]`` above the
        final break.  None may equal exactly -1.
    transform : `Transform`, optional
        Map between physical and sampling space, by default :class:`Identity`.
    """

    def __init__(self, breaks, alphas, transform=None):
        super().__init__(transform)
        self.breaks = np.asarray(breaks, dtype=float)
        self.alphas = np.asarray(alphas, dtype=float)
        if self.alphas.shape != (self.breaks.size + 1,):
            raise ValueError("alphas must have exactly one more entry than breaks.")
        if np.any(np.diff(self.breaks) <= 0):
            raise ValueError("breaks must be strictly increasing.")
        if np.any(self.alphas == -1.0):
            raise ValueError("BrokenPowerLaw does not support an exponent of exactly -1.")
        # Per-segment continuity coefficients (the lowest segment is 1).
        coeffs = np.ones(self.alphas.size)
        for i in range(self.breaks.size):
            coeffs[i + 1] = coeffs[i] * self.breaks[i] ** (self.alphas[i] - self.alphas[i + 1])
        self.coeffs = coeffs

    def _segments(self, lo, hi):
        """Decompose ``[lo, hi]`` into segments split at the in-range breaks.

        Returns the segment edges and, per segment, the exponent, continuity
        coefficient, and cumulative unnormalised integral (``cum[-1]`` is the
        normalisation constant ``Z``).
        """
        interior = self.breaks[(self.breaks > lo) & (self.breaks < hi)]
        edges = np.concatenate(([lo], interior, [hi]))
        piece = np.searchsorted(self.breaks, edges[:-1], side='right')
        alpha = self.alphas[piece]
        coeff = self.coeffs[piece]
        a = alpha + 1.0
        seg_int = coeff * (edges[1:] ** a - edges[:-1] ** a) / a
        cum = np.concatenate(([0.0], np.cumsum(seg_int)))
        return edges, alpha, coeff, cum

    def pdf(self, values, lo, hi):
        edges, alpha, coeff, cum = self._segments(lo, hi)
        s = np.clip(np.searchsorted(edges, values, side='right') - 1, 0, alpha.size - 1)
        return coeff[s] * np.power(values, alpha[s]) / cum[-1]

    def sample(self, n, lo, hi, rng=None):
        rng = rng or np.random.default_rng()
        return self._ppf(rng.uniform(0, 1, n), lo, hi)

    def sigma(self, values, lo, hi, avg_density):
        """CDF-space step, identical in spirit to :meth:`PowerLaw.sigma`."""
        u = self._cdf(values, lo, hi)
        x_right = self._ppf(np.clip(u + avg_density, 0.0, 1.0), lo, hi)
        x_left = self._ppf(np.clip(u - avg_density, 0.0, 1.0), lo, hi)
        return np.maximum(np.abs(x_right - values), np.abs(x_left - values))

    def _cdf(self, values, lo, hi):
        edges, alpha, coeff, cum = self._segments(lo, hi)
        s = np.clip(np.searchsorted(edges, values, side='right') - 1, 0, alpha.size - 1)
        a = alpha[s] + 1.0
        partial = coeff[s] * (np.power(values, a) - np.power(edges[s], a)) / a
        return (cum[s] + partial) / cum[-1]

    def _ppf(self, u, lo, hi):
        edges, alpha, coeff, cum = self._segments(lo, hi)
        target = np.clip(u, 0.0, 1.0) * cum[-1]
        s = np.clip(np.searchsorted(cum, target, side='right') - 1, 0, alpha.size - 1)
        a = alpha[s] + 1.0
        partial = target - cum[s]
        return np.power(partial * a / coeff[s] + np.power(edges[s], a), 1.0 / a)


class TruncatedNormal(Distribution):
    """Normal distribution truncated to ``[lo, hi]`` in sampling space.

    Combined with :class:`Ln` this gives the ``log_normal`` natal-kick prior:
    the kick magnitude ``v`` follows ``LogNormal(mu, scale)`` so ``ln(v)`` is
    normally distributed, and sampling space is ``ln(v)``.

    Parameters
    ----------
    mu : `float`
        Mean of the underlying normal (in sampling space).
    scale : `float`
        Standard deviation of the underlying normal (in sampling space).
    transform : `Transform`, optional
        Map between physical and sampling space, by default :class:`Identity`.
    """

    def __init__(self, mu, scale, transform=None):
        super().__init__(transform)
        self.mu = mu
        self.scale = scale

    def sample(self, n, lo, hi, rng=None):
        rng = rng or np.random.default_rng()
        # Truncated-normal inverse CDF: map uniform draws into
        # [CDF(lo), CDF(hi)] then apply the inverse normal CDF.
        p_lo = _scipy_norm.cdf(lo, loc=self.mu, scale=self.scale)
        p_hi = _scipy_norm.cdf(hi, loc=self.mu, scale=self.scale)
        u = rng.uniform(p_lo, p_hi, n)
        return _scipy_norm.ppf(u, loc=self.mu, scale=self.scale)

    def pdf(self, values, lo, hi):
        p_lo = _scipy_norm.cdf(lo, loc=self.mu, scale=self.scale)
        p_hi = _scipy_norm.cdf(hi, loc=self.mu, scale=self.scale)
        norm_factor = p_hi - p_lo  # probability mass within bounds
        return _scipy_norm.pdf(values, loc=self.mu, scale=self.scale) / norm_factor


# ---------------------------------------------------------------------------
# Registry
# ---------------------------------------------------------------------------
DISTRIBUTIONS = {
    'uniform': Uniform(),
    'flat_in_log': Uniform(transform=Log10()),
    'uniform_in_sine': Uniform(transform=Sin()),
    'uniform_in_cosine': Uniform(transform=CosShift()),
    'kroupa': BrokenPowerLaw(breaks=[0.5], alphas=[-1.3, -2.3]),
    'sana': PowerLaw(-0.55, transform=Log10()),
    'sana_ecc': PowerLaw(-0.45),
    'disberg': TruncatedNormal(5.67, 0.59, transform=Ln()),
}


def register(name, distribution):
    """Register a distribution instance under ``name``.

    Parameters
    ----------
    name : `str`
        Key used to refer to the distribution (e.g. from a
        :class:`~cosmic.sample.stroopwafel.parameter_space.Parameter`).
    distribution : `Distribution`
        The distribution instance to register.
    """
    if not isinstance(distribution, Distribution):
        raise TypeError(
            f"register() expects a Distribution instance, got "
            f"{type(distribution).__name__}."
        )
    DISTRIBUTIONS[name] = distribution


def get_distribution(dist):
    """Resolve a name or instance to a :class:`Distribution`.

    Parameters
    ----------
    dist : `str` or `Distribution`
        Either a key in :data:`DISTRIBUTIONS` or a distribution instance
        (returned unchanged).

    Returns
    -------
    `Distribution`
        The resolved distribution instance.
    """
    if isinstance(dist, Distribution):
        return dist
    try:
        return DISTRIBUTIONS[dist]
    except KeyError:
        raise KeyError(
            f"Unknown distribution {dist!r}.  Registered names: "
            f"{sorted(DISTRIBUTIONS)}.  Alternatively pass a Distribution "
            f"instance directly."
        ) from None
