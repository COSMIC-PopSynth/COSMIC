"""Vectorized rejection checking for binary star systems.

Provides fully vectorized numpy operations for ZAMS radius, Roche lobe
radius, and physical rejection criteria used to discard unphysical binary
systems prior to evolution.
"""
from cosmic.sample.sampler.independent import Sample
from cosmic.utils import calc_Roche_radius, a_from_p


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

    # compute separation from periods and masses
    separation = a_from_p(p=porb, m1=mass_1, m2=mass_2)

    # get stellar radii at ZAMS
    sampler = Sample()
    radius_1 = sampler.set_reff(mass=mass_1, metallicity=metallicity)
    radius_2 = sampler.set_reff(mass=mass_2, metallicity=metallicity)

    # roche lobe radii at periastron
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
