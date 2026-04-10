"""Vectorized prior probability density functions.

Each function takes ``(values, lo, hi)`` where ``values`` is an (N,) ndarray
and ``lo``/``hi`` are the bounds in the sampling (transformed) space, and
returns an (N,) ndarray of prior probability densities.
"""
import numpy as np
from .constants import ALPHA_IMF, SANA_G, SANA_ECC


def uniform(values, lo, hi):
    """Compute the uniform prior probability density.

    Parameters
    ----------
    values : `numpy.ndarray`
        (N,) array of sample values (unused, density is constant).
    lo : `float`
        Lower bound of the uniform distribution.
    hi : `float`
        Upper bound of the uniform distribution.

    Returns
    -------
    `numpy.ndarray`
        (N,) array of constant prior densities.
    """
    return np.full(len(values), 1.0 / (hi - lo))


def flat_in_log(values, lo, hi):
    """Compute the flat-in-log prior probability density.

    The prior is uniform in log10 space; bounds are already in log10.

    Parameters
    ----------
    values : `numpy.ndarray`
        (N,) array of sample values in log10 space (unused, density is
        constant).
    lo : `float`
        Lower bound in log10 space.
    hi : `float`
        Upper bound in log10 space.

    Returns
    -------
    `numpy.ndarray`
        (N,) array of constant prior densities.
    """
    return np.full(len(values), 1.0 / (hi - lo))


def kroupa(values, lo, hi):
    """Compute the Kroupa IMF power-law prior probability density.

    Parameters
    ----------
    values : `numpy.ndarray`
        (N,) array of sample values.
    lo : `float`
        Lower bound of the distribution.
    hi : `float`
        Upper bound of the distribution.

    Returns
    -------
    `numpy.ndarray`
        (N,) array of prior densities.
    """
    a = ALPHA_IMF
    norm = (a + 1) / (hi**(a + 1) - lo**(a + 1))
    return norm * np.power(values, a)


def sana(values, lo, hi):
    """Compute the Sana orbital period power-law prior probability density.

    Parameters
    ----------
    values : `numpy.ndarray`
        (N,) array of sample values.
    lo : `float`
        Lower bound of the distribution.
    hi : `float`
        Upper bound of the distribution.

    Returns
    -------
    `numpy.ndarray`
        (N,) array of prior densities.
    """
    a = SANA_G
    norm = (a + 1) / (hi**(a + 1) - lo**(a + 1))
    return norm * np.power(values, a)


def sana_ecc(values, lo, hi):
    """Compute the Sana eccentricity power-law prior probability density.

    Parameters
    ----------
    values : `numpy.ndarray`
        (N,) array of sample values.
    lo : `float`
        Lower bound of the distribution.
    hi : `float`
        Upper bound of the distribution.

    Returns
    -------
    `numpy.ndarray`
        (N,) array of prior densities.
    """
    a = SANA_ECC
    norm = (a + 1) / (hi**(a + 1) - lo**(a + 1))
    return norm * np.power(values, a)


def uniform_in_sine(values, lo, hi):
    """Compute the uniform-in-sine prior probability density.

    Parameters
    ----------
    values : `numpy.ndarray`
        (N,) array of sample values in sine space (unused, density is
        constant).
    lo : `float`
        Lower bound in sine space.
    hi : `float`
        Upper bound in sine space.

    Returns
    -------
    `numpy.ndarray`
        (N,) array of constant prior densities.
    """
    return np.full(len(values), 1.0 / (hi - lo))


def uniform_in_cosine(values, lo, hi):
    """Compute the uniform-in-cosine prior probability density.

    Parameters
    ----------
    values : `numpy.ndarray`
        (N,) array of sample values in cosine space (unused, density is
        constant).
    lo : `float`
        Lower bound in cosine space.
    hi : `float`
        Upper bound in cosine space.

    Returns
    -------
    `numpy.ndarray`
        (N,) array of constant prior densities.
    """
    return np.full(len(values), 1.0 / (hi - lo))


# Registry mapping string names to functions
PRIORS = {
    'uniform': uniform,
    'flat_in_log': flat_in_log,
    'kroupa': kroupa,
    'sana': sana,
    'sana_ecc': sana_ecc,
    'uniform_in_sine': uniform_in_sine,
    'uniform_in_cosine': uniform_in_cosine,
}
