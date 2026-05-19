"""Vectorized sampling functions for each distribution type.

Each function draws ``n`` samples between ``lo`` and ``hi`` bounds (in the
sampling/transformed space) and returns an (N,) ndarray.
"""
import numpy as np
from scipy.stats import norm as _scipy_norm
from .constants import ALPHA_IMF, SANA_G, SANA_ECC, NATAL_KICK_LOG_MU, NATAL_KICK_LOG_SIGMA


def uniform(n, lo, hi, rng=None):
    """Draw samples from a uniform distribution.

    Parameters
    ----------
    n : `int`
        Number of samples to draw.
    lo : `float`
        Lower bound of the uniform distribution.
    hi : `float`
        Upper bound of the uniform distribution.
    rng : `numpy.random.Generator`, optional
        Random number generator, by default None

    Returns
    -------
    `numpy.ndarray`
        (N,) array of uniform samples.
    """
    rng = rng or np.random.default_rng()
    return rng.uniform(lo, hi, n)


def flat_in_log(n, lo, hi, rng=None):
    """Sample uniformly in log10 space.

    Parameters
    ----------
    n : `int`
        Number of samples to draw.
    lo : `float`
        Lower bound (already in log10 space).
    hi : `float`
        Upper bound (already in log10 space).
    rng : `numpy.random.Generator`, optional
        Random number generator, by default None

    Returns
    -------
    `numpy.ndarray`
        (N,) array of samples in log10 space.
    """
    rng = rng or np.random.default_rng()
    return rng.uniform(lo, hi, n)


def kroupa(n, lo, hi, rng=None):
    """Inverse CDF sampling from a Kroupa-like power law p(x) ~ x^alpha.

    Parameters
    ----------
    n : `int`
        Number of samples to draw.
    lo : `float`
        Lower bound of the distribution.
    hi : `float`
        Upper bound of the distribution.
    rng : `numpy.random.Generator`, optional
        Random number generator, by default None

    Returns
    -------
    `numpy.ndarray`
        (N,) array of power-law distributed samples.
    """
    rng = rng or np.random.default_rng()
    u = rng.uniform(0, 1, n)
    a = ALPHA_IMF + 1
    return np.power(u * (hi**a - lo**a) + lo**a, 1.0 / a)


def sana(n, lo, hi, rng=None):
    """Inverse CDF sampling from the Sana orbital period distribution.

    The Sana et al. (2012) distribution is a power law in log10(period),
    so this sampler operates entirely in log10(period / days) space.
    ``lo`` and ``hi`` must therefore be given as log10 values (e.g.
    ``lo=0.15, hi=5.5`` spans ~1.4 d to ~316 000 d).

    Parameters
    ----------
    n : `int`
        Number of samples to draw.
    lo : `float`
        Lower bound in log10(period / days).
    hi : `float`
        Upper bound in log10(period / days).
    rng : `numpy.random.Generator`, optional
        Random number generator, by default None

    Returns
    -------
    `numpy.ndarray`
        (N,) array of samples from the Sana period distribution.
    """
    rng = rng or np.random.default_rng()
    u = rng.uniform(0, 1, n)
    a = SANA_G + 1
    return np.power(u * (hi**a - lo**a) + lo**a, 1.0 / a)


def sana_ecc(n, lo, hi, rng=None):
    """Inverse CDF sampling from the Sana eccentricity distribution.

    Parameters
    ----------
    n : `int`
        Number of samples to draw.
    lo : `float`
        Lower bound of the eccentricity distribution.
    hi : `float`
        Upper bound of the eccentricity distribution.
    rng : `numpy.random.Generator`, optional
        Random number generator, by default None

    Returns
    -------
    `numpy.ndarray`
        (N,) array of samples from the Sana eccentricity distribution.
    """
    rng = rng or np.random.default_rng()
    u = rng.uniform(0, 1, n)
    a = SANA_ECC + 1
    return np.power(u * (hi**a - lo**a) + lo**a, 1.0 / a)


def uniform_in_sine(n, lo, hi, rng=None):
    """Sample uniformly in sine-transformed space.

    Parameters
    ----------
    n : `int`
        Number of samples to draw.
    lo : `float`
        Lower bound in sine space.
    hi : `float`
        Upper bound in sine space.
    rng : `numpy.random.Generator`, optional
        Random number generator, by default None

    Returns
    -------
    `numpy.ndarray`
        (N,) array of uniform samples in sine space.
    """
    rng = rng or np.random.default_rng()
    return rng.uniform(lo, hi, n)


def uniform_in_cosine(n, lo, hi, rng=None):
    """Sample uniformly in cosine-transformed space.

    Parameters
    ----------
    n : `int`
        Number of samples to draw.
    lo : `float`
        Lower bound in cosine space.
    hi : `float`
        Upper bound in cosine space.
    rng : `numpy.random.Generator`, optional
        Random number generator, by default None

    Returns
    -------
    `numpy.ndarray`
        (N,) array of uniform samples in cosine space.
    """
    rng = rng or np.random.default_rng()
    return rng.uniform(lo, hi, n)


def log_normal(n, lo, hi, rng=None):
    """Inverse-CDF sampling from a (truncated) log-normal natal kick distribution.

    The kick magnitude v follows LogNormal(``NATAL_KICK_LOG_MU``,
    ``NATAL_KICK_LOG_SIGMA``), so ``ln(v)`` is normally distributed.
    ``lo`` and ``hi`` are bounds in ``ln(v)`` space (i.e. the natural log of
    the physical bounds in km/s), and sampling is restricted to that range
    via the truncated-normal inverse CDF.

    With ``NATAL_KICK_LOG_MU = 5.67`` and ``NATAL_KICK_LOG_SIGMA = 0.59``,
    the median kick is ``exp(5.67) ≈ 291 km/s`` (Hobbs et al. 2005 / Fryer
    et al.).  Physical bounds of ``[0.1, 5000] km/s`` map to
    ``lo ≈ –2.30``, ``hi ≈ 8.52`` in sampling space, which captures
    essentially all of the probability mass.

    Parameters
    ----------
    n : `int`
        Number of samples to draw.
    lo : `float`
        Lower bound in ``ln(v / km s⁻¹)`` space.
    hi : `float`
        Upper bound in ``ln(v / km s⁻¹)`` space.
    rng : `numpy.random.Generator`, optional
        Random number generator, by default None

    Returns
    -------
    `numpy.ndarray`
        (N,) array of samples in ``ln(v)`` space.
    """
    rng = rng or np.random.default_rng()
    # Truncated-normal inverse CDF: map uniform draws to [CDF(lo), CDF(hi)],
    # then apply the inverse normal CDF.
    p_lo = _scipy_norm.cdf(lo, loc=NATAL_KICK_LOG_MU, scale=NATAL_KICK_LOG_SIGMA)
    p_hi = _scipy_norm.cdf(hi, loc=NATAL_KICK_LOG_MU, scale=NATAL_KICK_LOG_SIGMA)
    u = rng.uniform(p_lo, p_hi, n)
    return _scipy_norm.ppf(u, loc=NATAL_KICK_LOG_MU, scale=NATAL_KICK_LOG_SIGMA)


# Registry mapping string names to functions
SAMPLERS = {
    'uniform': uniform,
    'flat_in_log': flat_in_log,
    'kroupa': kroupa,
    'sana': sana,
    'sana_ecc': sana_ecc,
    'uniform_in_sine': uniform_in_sine,
    'uniform_in_cosine': uniform_in_cosine,
    'log_normal': log_normal,
}
