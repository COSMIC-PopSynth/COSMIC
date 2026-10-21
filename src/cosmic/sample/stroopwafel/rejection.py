"""Vectorized rejection checking for binary star systems.

Provides fully vectorized numpy operations for ZAMS radius, Roche lobe
radius, and physical rejection criteria used to discard unphysical binary
systems prior to evolution.
"""
from functools import lru_cache

import numpy as np

from cosmic.sample.sampler.independent import Sample
from cosmic.utils import calc_Roche_radius, a_from_p


def _using_metisse(SSEDict):
    return SSEDict is not None and SSEDict.get("stellar_engine") == "metisse"


@lru_cache(maxsize=8)
def _metisse_zams_grid(metallicity, path_to_tracks, path_to_he_tracks, z_accuracy_limit, n_grid):
    """Tabulate METISSE ZAMS radii on a log-spaced grid spanning the track masses

    Cached since the grid only depends on the tracks and metallicity, so it is computed once per
    process rather than for every call to the rejection function

    Returns
    -------
    log_mass, log_radius : `numpy.ndarray`
        (n_grid,) arrays of the natural log of the grid masses and their ZAMS radii
    """
    from cosmic.evolve import read_tracks_for_METISSE

    m_min, m_max = read_tracks_for_METISSE(path_to_tracks=path_to_tracks, IBT_Z=metallicity,
                                           z_accuracy_limit=z_accuracy_limit, is_he=False)
    SSEDict = {"stellar_engine": "metisse", "path_to_tracks": path_to_tracks,
               "path_to_he_tracks": path_to_he_tracks, "z_accuracy_limit": z_accuracy_limit}
    grid = np.geomspace(m_min, m_max, n_grid)
    radii = Sample().set_reff(mass=grid, metallicity=metallicity, SSEDict=SSEDict)
    return np.log(grid), np.log(radii)


def zams_radii(mass, metallicity, SSEDict=None, n_grid=5000):
    """ZAMS radii for a set of masses at a single metallicity

    With METISSE, computing a ZAMS radius costs ~0.1 ms per star (interpolating a full track), which
    dominates rejection when it is applied to millions of samples.  Since the radius depends only on
    the mass (for a given metallicity and set of tracks), METISSE radii are instead tabulated once on
    a log-spaced grid of ``n_grid`` masses spanning the hydrogen tracks and then linearly interpolated
    in log-log space (relative errors ~1e-4 for the default grid).  Masses outside the grid's range
    are computed directly.  SSE radii are cheap and are always computed directly.

    Parameters
    ----------
    mass : `numpy.ndarray`
        (N,) array of masses in solar masses
    metallicity : `float`
        Metallicity of all stars
    SSEDict : `dict`, optional
        COSMIC single stellar evolution settings, by default None (COSMIC's ``sse`` engine)
    n_grid : `int` or None, optional
        Number of grid masses used for METISSE radii, by default 5000.  Pass None to compute every
        radius directly with ``Sample.set_reff``

    Returns
    -------
    `numpy.ndarray`
        (N,) array of ZAMS radii in solar radii
    """
    mass = np.asarray(mass, dtype=float)
    if n_grid is None or not _using_metisse(SSEDict):
        return Sample().set_reff(mass=mass, metallicity=metallicity, SSEDict=SSEDict)

    log_m, log_r = _metisse_zams_grid(float(metallicity), SSEDict["path_to_tracks"],
                                      SSEDict["path_to_he_tracks"],
                                      SSEDict.get("z_accuracy_limit", 1e-2), int(n_grid))
    log_mass = np.log(mass)
    radii = np.exp(np.interp(log_mass, log_m, log_r))

    outside = (log_mass < log_m[0]) | (log_mass > log_m[-1])
    if np.any(outside):
        radii[outside] = Sample().set_reff(mass=mass[outside], metallicity=metallicity,
                                           SSEDict=SSEDict)
    return radii


def default_reject(binary_params, SSEDict=None, min_secondary_mass=0.08, n_radius_grid=5000):
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
    SSEDict : `dict`, optional
        COSMIC single stellar evolution settings.  The ZAMS radii are
        obtained from ``Sample.set_reff``, which depends on the stellar
        engine (e.g. ``sse`` vs ``METISSE``); passing the same ``SSEDict``
        used for evolution ensures the rejection radii are computed
        consistently.  By default None (COSMIC's ``sse`` engine).
        :class:`~cosmic.sample.stroopwafel.main.AdaptiveSampler` wires its
        own ``SSEDict`` into this argument automatically.
    min_secondary_mass : `float`, optional
        Minimum allowed secondary mass in solar masses, by default 0.08
    n_radius_grid : `int` or None, optional
        Number of grid masses used to tabulate METISSE ZAMS radii, see
        :func:`zams_radii`.  Pass None to compute every radius directly.
        By default 5000

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

    # if using METISSE, the metallicities must be a single value (not an array) to compute ZAMS radii
    if _using_metisse(SSEDict):
        if not isinstance(metallicity, float):
            if len(set(metallicity)) > 1:
                raise ValueError(
                    "METISSE requires a single metallicity value to compute ZAMS radii, "
                    "but multiple metallicities were provided."
                )
            metallicity = metallicity[0]

    # get stellar radii at ZAMS (depends on the stellar engine in SSEDict),
    # using a single call since each call (re)loads the stellar tracks
    radii = zams_radii(np.concatenate([mass_1, mass_2]), metallicity,
                       SSEDict=SSEDict, n_grid=n_radius_grid)
    radius_1, radius_2 = radii[:len(mass_1)], radii[len(mass_1):]

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
