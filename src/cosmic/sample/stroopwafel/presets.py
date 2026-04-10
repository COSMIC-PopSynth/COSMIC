"""Preset hit-definition functions for common DCO selections.

Each preset factory returns a callable with signature
``(bpp) -> (n_hits, hit_bin_nums)`` where ``hit_bin_nums`` are the
``bin_num`` values of systems classified as hits.
"""
import numpy as np


def merging_dco(kstar_1, kstar_2, max_merge_time=13.7):
    """Create a hit function selecting merging double compact objects.

    Parameters
    ----------
    kstar_1 : `list` of `int`
        Allowed kstar types for star 1 (e.g., ``[14]`` for black holes).
    kstar_2 : `list` of `int`
        Allowed kstar types for star 2 (e.g., ``[14]`` for black holes).
    max_merge_time : `float`, optional
        Maximum merger time in Gyr, by default 13.7 (Hubble time)

    Returns
    -------
    `callable`
        Function with signature ``(bpp) -> (n_hits, hit_bin_nums)``
        where ``bpp`` is a `pandas.DataFrame` and ``hit_bin_nums`` is a
        `numpy.ndarray` of bin_num values.
    """
    from legwork import evol
    import astropy.units as u

    k1_set = set(kstar_1)
    k2_set = set(kstar_2)

    def is_interesting(bpp):
        # Select rows matching the DCO type (either ordering)
        pairs_mask = (
            (bpp.kstar_1.isin(k1_set) & bpp.kstar_2.isin(k2_set))
            | (bpp.kstar_1.isin(k2_set) & bpp.kstar_2.isin(k1_set))
        )
        # Must still be bound (sep > 0)
        interesting_mask = pairs_mask & (bpp.sep > 0)
        candidates = bpp.loc[interesting_mask].drop_duplicates(subset='bin_num', keep='first')

        if len(candidates) == 0:
            return 0, np.array([], dtype=int)

        # Compute merger times using LEGWORK
        merge_times = evol.get_t_merge_ecc(
            ecc_i=candidates.ecc.values,
            a_i=candidates.sep.values * u.Rsun,
            m_1=candidates.mass_1.values * u.Msun,
            m_2=candidates.mass_2.values * u.Msun,
        )
        merge_mask = merge_times < (max_merge_time * u.Gyr)
        hits = candidates[merge_mask]

        return len(hits), hits.bin_num.values

    return is_interesting


def any_dco(kstar_1, kstar_2):
    """Create a hit function selecting DCOs regardless of merge time.

    Parameters
    ----------
    kstar_1 : `list` of `int`
        Allowed kstar types for star 1.
    kstar_2 : `list` of `int`
        Allowed kstar types for star 2.

    Returns
    -------
    `callable`
        Function with signature ``(bpp) -> (n_hits, hit_bin_nums)``
        where ``bpp`` is a `pandas.DataFrame` and ``hit_bin_nums`` is a
        `numpy.ndarray` of bin_num values.
    """
    k1_set = set(kstar_1)
    k2_set = set(kstar_2)

    def is_interesting(bpp):
        pairs_mask = (
            (bpp.kstar_1.isin(k1_set) & bpp.kstar_2.isin(k2_set))
            | (bpp.kstar_1.isin(k2_set) & bpp.kstar_2.isin(k1_set))
        )
        interesting_mask = pairs_mask & (bpp.sep > 0)
        candidates = bpp.loc[interesting_mask].drop_duplicates(subset='bin_num', keep='first')
        return len(candidates), candidates.bin_num.values

    return is_interesting
