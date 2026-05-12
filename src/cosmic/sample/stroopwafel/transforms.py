"""Vectorized scale transformations between physical and sampling space.

Some distributions (flat_in_log, sana, uniform_in_sine, etc.) sample in a
transformed space (e.g., log10). These functions convert between the two
representations. All functions operate on (N,) column arrays for a given
parameter.
"""
import numpy as np


def to_sampling_space(values, sampler_type):
    """Convert physical-space values to the sampling (transformed) space.

    Parameters
    ----------
    values : `numpy.ndarray`
        (N,) array of values in physical space.
    sampler_type : `str`
        Name of the sampling distribution (e.g., ``'flat_in_log'``,
        ``'uniform_in_sine'``).

    Returns
    -------
    `numpy.ndarray`
        (N,) array of values in sampling space.
    """
    if sampler_type in ('flat_in_log', 'sana'):
        return np.log10(values)
    elif sampler_type == 'uniform_in_sine':
        return np.sin(values)
    elif sampler_type == 'uniform_in_cosine':
        # Angles are measured from –π/2 to π/2 (e.g. declination-like
        # coordinates), so the sampling variable is cos(θ + π/2) = –sin(θ).
        # The round-trip is exact for θ ∈ [–π/2, π/2].
        return np.cos(values + np.pi / 2)
    return values


def to_physical_space(values, sampler_type):
    """Convert sampling-space values back to physical space.

    Parameters
    ----------
    values : `numpy.ndarray`
        (N,) array of values in sampling space.
    sampler_type : `str`
        Name of the sampling distribution (e.g., ``'flat_in_log'``,
        ``'uniform_in_sine'``).

    Returns
    -------
    `numpy.ndarray`
        (N,) array of values in physical space.
    """
    if sampler_type in ('flat_in_log', 'sana'):
        return np.power(10.0, values)
    elif sampler_type == 'uniform_in_sine':
        return np.arcsin(values)
    elif sampler_type == 'uniform_in_cosine':
        # Inverse of cos(θ + π/2): arccos(u) – π/2
        return np.arccos(values) - np.pi / 2
    return values


def transform_bounds(lo, hi, sampler_type):
    """Get the bounds in sampling space for a given physical-space range.

    Parameters
    ----------
    lo : `float`
        Lower bound in physical space.
    hi : `float`
        Upper bound in physical space.
    sampler_type : `str`
        Name of the sampling distribution.

    Returns
    -------
    lo_transformed : `float`
        Lower bound in sampling space.
    hi_transformed : `float`
        Upper bound in sampling space.
    """
    if sampler_type == 'flat_in_log':
        return np.log10(lo), np.log10(hi)
    elif sampler_type == 'uniform_in_sine':
        return -1.0, 1.0
    elif sampler_type == 'uniform_in_cosine':
        return -1.0, 1.0
    # kroupa, uniform: bounds are already in physical space.
    # sana, sana_ecc: bounds are passed by the caller in the native
    # sampling space (log10(period) for sana, eccentricity for sana_ecc),
    # so no further transformation is needed here.
    return lo, hi
