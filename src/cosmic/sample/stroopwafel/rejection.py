"""Vectorized rejection checking for binary star systems.

Provides fully vectorized numpy operations for ZAMS radius, Roche lobe
radius, and physical rejection criteria used to discard unphysical binary
systems prior to evolution.
"""
import numpy as np
from .constants import R_COEFF, ZSOL, R_SOL_TO_AU
from cosmic.utils import calc_Roche_radius


def get_zams_radius(mass, metallicity):
    """Compute zero-age main sequence radius for arrays of stars.

    Parameters
    ----------
    mass : `numpy.ndarray`
        (N,) array of stellar masses in solar masses.
    metallicity : `numpy.ndarray`
        (N,) array of metallicities.

    Returns
    -------
    `numpy.ndarray`
        (N,) array of ZAMS radii in AU.
    """
    mass = np.asarray(mass, dtype=float)
    metallicity = np.asarray(metallicity, dtype=float)

    xi = np.log10(metallicity / ZSOL)

    # R_COEFF is a (9, 5) matrix of polynomial coefficients
    # For each of the 9 radius coefficients, evaluate the polynomial in xi
    R = np.array(R_COEFF)  # (9, 5)
    # Build Vandermonde matrix: xi^0, xi^1, xi^2, xi^3, xi^4
    powers = np.column_stack([xi**k for k in range(5)])  # (N, 5)
    # rc[j, n] = sum_k R[j, k] * xi[n]^k
    rc = R @ powers.T  # (9, N)

    top = (rc[0] * mass**2.5 + rc[1] * mass**6.5 + rc[2] * mass**11
           + rc[3] * mass**19 + rc[4] * mass**19.5)
    bottom = (rc[5] + rc[6] * mass**2 + rc[7] * mass**8.5
              + mass**18.5 + rc[8] * mass**19.5)

    return (top / bottom) * R_SOL_TO_AU


def default_reject(binary_params, min_secondary_mass=0.08):
    """Default rejection function for DCO progenitor systems.

    Rejects systems where the secondary mass is below the minimum, the
    stars are in contact at ZAMS, or either star overflows its Roche lobe
    at periastron.  The orbital separation is computed from the orbital
    period via Kepler's third law, and both stars share the binary
    metallicity.

    Parameters
    ----------
    binary_params : `dict`
        Assembled binary parameters with keys ``'mass_1'``, ``'mass_2'``
        (solar masses), ``'porb'`` (days), ``'ecc'``, and ``'metallicity'``,
        each an (N,) array.
    min_secondary_mass : `float`, optional
        Minimum allowed secondary mass in solar masses, by default 0.08

    Returns
    -------
    `numpy.ndarray`
        (N,) boolean mask where True indicates a rejected system.
    """
    mass_1 = binary_params['mass_1']
    mass_2 = binary_params['mass_2']
    porb = binary_params['porb']
    ecc = binary_params['ecc']
    metallicity = binary_params['metallicity']

    # Semi-major axis [AU] from Kepler's third law (porb in days, masses in
    # solar masses): a^3 [AU^3] = (P [yr])^2 * M [Msun].
    separation = ((porb / 365.25) ** 2 * (mass_1 + mass_2)) ** (1.0 / 3.0)

    # ZAMS radii [AU] (both stars share the binary metallicity)
    radius_1 = get_zams_radius(mass_1, metallicity)
    radius_2 = get_zams_radius(mass_2, metallicity)

    # Roche lobe radii at periastron
    peri_sep = separation * (1 - ecc)
    rl_1 = calc_Roche_radius(mass_1, mass_2, peri_sep)
    rl_2 = calc_Roche_radius(mass_2, mass_1, peri_sep)

    rejected = (
        (mass_2 < min_secondary_mass)
        | (separation <= (radius_1 + radius_2))
        | (radius_1 / rl_1 > 1)
        | (radius_2 / rl_2 > 1)
    )

    return rejected
