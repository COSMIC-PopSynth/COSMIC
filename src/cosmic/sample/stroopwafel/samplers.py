"""Vectorized sampling functions for each distribution type.

Each function draws ``n`` samples between ``lo`` and ``hi`` bounds (in the
sampling/transformed space) and returns an (N,) ndarray.
"""
import numpy as np
from .constants import ALPHA_IMF, SANA_G, SANA_ECC


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

    Parameters
    ----------
    n : `int`
        Number of samples to draw.
    lo : `float`
        Lower bound (in log10 space).
    hi : `float`
        Upper bound (in log10 space).
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


# Registry mapping string names to functions
SAMPLERS = {
    'uniform': uniform,
    'flat_in_log': flat_in_log,
    'kroupa': kroupa,
    'sana': sana,
    'sana_ecc': sana_ecc,
    'uniform_in_sine': uniform_in_sine,
    'uniform_in_cosine': uniform_in_cosine,
}
