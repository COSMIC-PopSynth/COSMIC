# -*- coding: utf-8 -*-
# Copyright (C) Katelyn Breivik (2017 - 2021)
#
# This file is part of cosmic.
#
# cosmic is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# cosmic is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with cosmic.  If not, see <http://www.gnu.org/licenses/>.

"""`independent`
"""

import numpy as np

from ... import utils

from .sampler import register_sampler
from .. import InitialBinaryTable
from ... import _evolvebin


__author__ = "Katelyn Breivik <katie.breivik@gmail.com>"
__credits__ = "Scott Coughlin <scott.coughlin@ligo.org>"
__all__ = ["get_independent_sampler", "Sample"]


def get_independent_sampler(
    final_kstar1,
    final_kstar2,
    primary_model,
    ecc_model,
    porb_model,
    SF_start,
    SF_duration,
    binfrac_model,
    met,
    size,
    **kwargs
):
    """Generates an initial binary sample according to user specified models

    Parameters
    ----------
    final_kstar1 : `int or list`
        Int or list of final kstar1

    final_kstar2 : `int or list`
        Int or list of final kstar2

    primary_model : `str`
        Model to sample primary mass; choices include: kroupa93, kroupa01, salpeter55, custom
        if 'custom' is selected, must also pass arguemts:
        alphas : `array`
            list of power law indicies
        mcuts : `array`
            breaks in the power laws.
        e.g. alphas=[-1.3,-2.3,-2.3],mcuts=[0.08,0.5,1.0,150.] reproduces standard Kroupa2001 IMF

    ecc_model : `str`
        Model to sample eccentricity; choices include: thermal, uniform, sana12

    porb_model : `str`
        Model to sample orbital period; choices include: log_uniform, sana12

    msort : `float`
        Stars with M>msort can have different pairing and sampling of companions

    pair : `float`
        Sets the pairing of stars M>msort only with stars with M>msort

    qmin : `float`
        kwarg which sets the minimum mass ratio for sampling the secondary
        where the mass ratio distribution is flat in q
        if q > 0, qmin sets the minimum mass ratio
        q = -1, this limits the minimum mass ratio to be set such that
        the pre-MS lifetime of the secondary is not longer than the full
        lifetime of the primary if it were to evolve as a single star

    m2_min : `float`
        kwarg which sets the minimum secondary mass for sampling
        the secondary as uniform in mass_2 between m2_min and mass_1

    qmin_msort : `float`
        Same as qmin for M>msort; only applies if qmin is supplied

    SF_start : `float`
        Time in the past when star formation initiates in Myr

    SF_duration : `float`
        Duration of constant star formation beginning from SF_Start in Myr

    binfrac_model : `str or float`
        Model for binary fraction; choices include: vanHaaften or a fraction where 1.0 is 100% binaries

    binfrac_model_msort : `str or float`
        Same as binfrac_model for M>msort

    met : `float`
        Sets the metallicity of the binary population where solar metallicity is zsun 

    size : `int`
        Size of the population to sample

    zsun : `float`
        optional kwarg for setting effective radii, default is 0.02

    Returns
    -------
    InitialBinaryTable : `pandas.DataFrame`
        DataFrame in the format of the InitialBinaryTable

    mass_singles : `float`
        Total mass in single stars needed to generate population

    mass_binaries : `float`
        Total mass in binaries needed to generate population

    n_singles : `int`
        Number of single stars needed to generate a population

    n_binaries : `int`
        Number of binaries needed to generate a population
    """
    if type(final_kstar1) in [int, float]:
        final_kstar1 = [final_kstar1]
    if type(final_kstar2) in [int, float]:
        final_kstar2 = [final_kstar2]
    primary_min, primary_max, secondary_min, secondary_max = utils.mass_min_max_select(
        final_kstar1, final_kstar2
    )
    initconditions = Sample()

    # set up multiplier if the mass sampling is inefficient
    multiplier = 1
    mass1_singles = []
    mass1_binary = []
    mass2_binary = []
    binfrac = []

    # track the mass in singles and the mass in binaries
    mass_singles = 0.0
    mass_binaries = 0.0

    # track the total number of stars sampled
    n_singles = 0
    n_binaries = 0

    while len(mass1_binary) < size:
        mass1, total_mass1 = initconditions.sample_primary(
            primary_model, size=size * multiplier, **kwargs 
        )
        (
            mass1_binaries,
            mass_single,
            binfrac_binaries,
            binary_index,
        ) = initconditions.binary_select(mass1, binfrac_model=binfrac_model, **kwargs)
        mass2_binaries = initconditions.sample_secondary(
            mass1_binaries, **kwargs)

        # track the mass sampled
        mass_singles += np.sum(mass_single)
        mass_binaries += np.sum(mass1_binaries)
        mass_binaries += np.sum(mass2_binaries)

        # track the total number sampled
        n_singles += len(mass_single)
        n_binaries += len(mass1_binaries)

        # select out the primaries and secondaries that will produce the final kstars
        (ind_select_primary,) = np.where(
            (mass1_binaries > primary_min) & (mass1_binaries < primary_max)
        )
        (ind_select_secondary,) = np.where(
            (mass2_binaries > secondary_min) & (mass2_binaries < secondary_max)
        )
        ind_select = list(
            set(ind_select_primary).intersection(ind_select_secondary))
        mass1_binary.extend(mass1_binaries[ind_select])
        mass2_binary.extend(mass2_binaries[ind_select])
        binfrac.extend(binfrac_binaries[ind_select])

        # select out the single stars that will produce the final kstar
        (ind_select_single,) = np.where(
            (mass_single > primary_min) & (mass_single < primary_max)
        )
        mass1_singles.extend(mass_single[ind_select_single])

        # check to see if we should increase the multiplier factor to sample the population more quickly
        if len(mass1_binary) < size / 100:
            # well this size clearly is not working time to increase
            # the multiplier by an order of magnitude
            multiplier *= 10

    mass1_binary = np.array(mass1_binary)
    mass2_binary = np.array(mass2_binary)
    binfrac = np.asarray(binfrac)
    mass1_singles = np.asarray(mass1_singles)

    zsun = kwargs.pop("zsun", 0.02)

    rad1 = initconditions.set_reff(mass1_binary, metallicity=met, zsun=zsun)
    rad2 = initconditions.set_reff(mass2_binary, metallicity=met, zsun=zsun)

    # sample periods and eccentricities
    porb,aRL_over_a = initconditions.sample_porb(
        mass1_binary, mass2_binary, rad1, rad2, porb_model, size=mass1_binary.size
    )
    ecc = initconditions.sample_ecc(aRL_over_a, ecc_model, size=mass1_binary.size)

    tphysf, metallicity = initconditions.sample_SFH(
        SF_start=SF_start, SF_duration=SF_duration, met=met, size=mass1_binary.size
    )
    metallicity[metallicity < 1e-4] = 1e-4
    metallicity[metallicity > 0.03] = 0.03
    kstar1 = initconditions.set_kstar(mass1_binary)
    kstar2 = initconditions.set_kstar(mass2_binary)

    if kwargs.pop("keep_singles", False):
        binary_table = InitialBinaryTable.InitialBinaries(
            mass1_binary,
            mass2_binary,
            porb,
            ecc,
            tphysf,
            kstar1,
            kstar2,
            metallicity,
            binfrac=binfrac,
        )
        tphysf, metallicity = initconditions.sample_SFH(
            SF_start=SF_start, SF_duration=SF_duration, met=met, size=mass1_singles.size
        )
        metallicity[metallicity < 1e-4] = 1e-4
        metallicity[metallicity > 0.03] = 0.03
        kstar1 = initconditions.set_kstar(mass1_singles)
        singles_table = InitialBinaryTable.InitialBinaries(
            mass1_singles,
            np.ones_like(mass1_singles)*0,
            np.ones_like(mass1_singles)*-1,
            np.ones_like(mass1_singles)*-1,
            tphysf,
            kstar1,
            np.ones_like(mass1_singles)*0,
            metallicity,
        )
        binary_table = binary_table.append(singles_table)
    else:
        binary_table = InitialBinaryTable.InitialBinaries(
            mass1_binary,
            mass2_binary,
            porb,
            ecc,
            tphysf,
            kstar1,
            kstar2,
            metallicity,
            binfrac=binfrac,
        )

    return (
        binary_table,
        mass_singles,
        mass_binaries,
        n_singles,
        n_binaries,
    )


register_sampler(
    "independent",
    InitialBinaryTable,
    get_independent_sampler,
    usage="final_kstar1, final_kstar2, binfrac_model, primary_model, ecc_model, SFH_model, component_age, metallicity, size",
)


class Sample(object):
    # sample primary masses
    def sample_primary(self, primary_model='kroupa01', size=None, **kwargs):
        """Sample the primary mass (always the most massive star) from a user-selected model

            kroupa93 follows Kroupa (1993), normalization comes from
            `Hurley 2002 <https://arxiv.org/abs/astro-ph/0201220>`_
            between 0.08 and 150 Msun
            salpter55 follows
            `Salpeter (1955) <http://adsabs.harvard.edu/abs/1955ApJ...121..161S>`_
            between 0.08 and 150 Msun
            kroupa01 follows Kroupa (2001) <https://arxiv.org/abs/astro-ph/0009005>
            between 0.08 and 100 Msun


            Parameters
            ----------
            primary_model : str, optional
            model for mass distribution; choose from:

            kroupa93 follows Kroupa (1993), normalization comes from
            `Hurley 2002 <https://arxiv.org/abs/astro-ph/0201220>`_
            valid for masses between 0.1 and 100 Msun

            salpter55 follows
            `Salpeter (1955) <http://adsabs.harvard.edu/abs/1955ApJ...121..161S>`_
            valid for masses between 0.1 and 100 Msun

            kroupa01 follows Kroupa (2001), normalization comes from
            `Hurley 2002 <https://arxiv.org/abs/astro-ph/0009005>`_
            valid for masses between 0.1 and 100 Msun

            custom is a generic piecewise power law that takes in the power
            law slopes and break points given in the optional input lists (alphas, mcuts)
            default alphas and mcuts yield an IMF identical to kroupa01

            Default kroupa01

            size : int, optional
            number of initial primary masses to sample
            NOTE: this is set in cosmic-pop call as Nstep

            alphas : array, optional
            absolute values of the power law slopes for primary_model = 'custom'
            Default [-1.3,-2.3,-2.3] (identical to slopes for primary_model = 'kroupa01')

            mcuts : array, optional, units of Msun
            break points separating the power law 'pieces' for primary_model = 'custom'
            Default [0.08,0.5,1.0,150.] (identical to breaks for primary_model = 'kroupa01')

            Returns
            -------
            a_0 : array
            Sampled primary masses
            np.sum(a_0) : float
            Total amount of mass sampled
            """


        if primary_model == 'kroupa93': 
            alphas, mcuts = [-1.3,-2.2,-2.7], [0.08,0.5,1.0,150.]
        # Since COSMIC/BSE can't handle < 0.08Msun, we will truncate at 0.08 Msun instead of 0.01
        elif primary_model == 'kroupa01': 
            alphas, mcuts = [-1.3,-2.3], [0.08,0.5,150.]
        elif primary_model == 'salpeter55': 
            alphas, mcuts = [-2.35], [0.08,150.]
        elif primary_model == 'custom': 
            if 'alphas' in kwargs and 'mcuts' in kwargs:
                alphas = kwargs.pop("alphas", [-1.3,-2.3,-2.3])
                mcuts = kwargs.pop("mcuts", [0.08,0.5,1.0,150.])
            else:
                raise ValueError("You must supply both alphas and mcuts to use"
                                 " a custom IMF generator")

        Ncumulative, Ntotal, coeff = [], 0., 1.
        for i in range(len(alphas)):
            g = 1. + alphas[i]
            # Compute this piece of the IMF's contribution to Ntotal
            if alphas[i] == -1: Ntotal += coeff * np.log(mcuts[i+1]/mcuts[i])
            else: Ntotal += coeff/g * (mcuts[i+1]**g - mcuts[i]**g)
            Ncumulative.append(Ntotal)
            if i < len(alphas)-1: coeff *= mcuts[i+1]**(-alphas[i+1]+alphas[i])

        cutoffs = np.array(Ncumulative)/Ntotal
        u = np.random.uniform(0.,1.,size)
        idxs = [() for i in range(len(alphas))]

        for i in range(len(alphas)):
            if i == 0: idxs[i], = np.where(u <= cutoffs[0])
            elif i < len(alphas)-1: idxs[i], = np.where((u > cutoffs[i-1]) & (u <= cutoffs[i]))
            else: idxs[i], = np.where(u > cutoffs[i-1])
        for i in range(len(alphas)):
            if alphas[i] == -1.0:
                u[idxs[i]] = 10**np.random.uniform(np.log10(mcuts[i]), 
                                                   np.log10(mcuts[i+1]), 
                                                   len(idxs[i]))
            else:
                u[idxs[i]] = utils.rndm(a=mcuts[i], b=mcuts[i+1], g=alphas[i], size=len(idxs[i]))

        return u, np.sum(u)

    # sample secondary mass
    def sample_secondary(self, primary_mass, **kwargs):
        """Sample a secondary mass using draws from a uniform mass ratio distribution motivated by
        `Mazeh et al. (1992) <http://adsabs.harvard.edu/abs/1992ApJ...401..265M>`_
        and `Goldberg & Mazeh (1994) <http://adsabs.harvard.edu/abs/1994ApJ...429..362G>`_

        NOTE: the lower lim is set by either qmin or m2_min which are passed as kwargs

        Parameters
        ----------
        primary_mass : `array`
            sets the maximum secondary mass (for a maximum mass ratio of 1)

        Returns
        -------
        secondary_mass : array
            sampled secondary masses with array size matching size of
            primary_mass
        """

        flag_msort = kwargs.pop("flag_msort", "no_msort")
        try:
            qmin = kwargs.pop("qmin")
        except:
            qmin = None
        try:
            m2_min = kwargs.pop("m2_min")
        except:
            m2_min = None

        if flag_msort == "binfrac_high_mass":
            msort = kwargs.pop("msort")
            qmin_msort = kwargs.pop("qmin_msort")
            pair = kwargs.pop("pair")
        else:
            msort, qmin_msort, pair = 10000., 0., 0

        sorting = msort * np.ones(primary_mass.size)
        (mhighIdx,) = np.where(primary_mass >= sorting)
        (mlowIdx,) = np.where(primary_mass < sorting)
        primary_mass_high = primary_mass[mhighIdx]
        primary_mass_low = primary_mass[mlowIdx]

        
        if qmin is not None:
            if (qmin >= 0.0):
                secondary_mass_low = np.random.uniform(
                    qmin, 1.0, len(primary_mass_low)
                ) * primary_mass_low
            elif (qmin < 0.0):
                dat = np.array([[5.0, 0.1363522012578616],
                                [6.999999999999993, 0.1363522012578616],
                                [12.599999999999994, 0.11874213836477984],
                                [20.999999999999993, 0.09962264150943395],
                                [29.39999999999999, 0.0820125786163522],
                                [41, 0.06490566037735851],
                                [55, 0.052327044025157254],
                                [70.19999999999999, 0.04301886792452836],
                                [87.4, 0.03622641509433966],
                                [107.40000000000002, 0.030188679245283068],
                                [133.40000000000003, 0.02515723270440262],
                                [156.60000000000002, 0.02163522012578628],
                                [175.40000000000003, 0.01962264150943399],
                                [200.20000000000005, 0.017358490566037776]])
                from scipy.interpolate import interp1d
                qmin_interp = interp1d(dat[:, 0], dat[:, 1])
                qmin = np.ones_like(primary_mass_low) * 0.1
                ind_5, = np.where(primary_mass_low > 5.0)
                qmin[ind_5] = qmin_interp(primary_mass_low[ind_5])

                q = np.random.uniform(qmin, 1.0)
                secondary_mass_low = q * primary_mass_low

        if m2_min is not None:
            m2_min = np.ones_like(primary_mass_low) * m2_min
            ind_m1_lim, = np.where(m2_min > primary_mass_low)
            m2_min[ind_m1_lim] = primary_mass_low[ind_m1_lim]
            secondary_mass_low = np.random.uniform(m2_min, primary_mass_low)

        if (m2_min is None) & (qmin is None):
            raise ValueError("You must supply either qmin or m2_min to the"
                             " independent initial binary sampler")

        if (msort < 10000.0) & (m2_min is None):
            if qmin_msort > 0.0 and pair == 0:
                secondary_mass_high = np.random.uniform(
                    qmin_msort, 1.0, len(primary_mass_high)
                ) * primary_mass_high
            elif qmin_msort > 0.0 and pair == 1:
                qmin_m = msort/primary_mass_high
                qmin_msort_arr = qmin_msort * np.ones(primary_mass_high.size)
                (qIdx,) = np.where(qmin_m < qmin_msort_arr)
                qmin_m[qIdx] = qmin_msort * np.ones(qIdx.size)
                q_msort = np.random.uniform(qmin_m, 1.0)
                secondary_mass_high = q_msort * primary_mass_high
            else:
                dat = np.array([[5.0, 0.1363522012578616],
                                [6.999999999999993, 0.1363522012578616],
                                [12.599999999999994, 0.11874213836477984],
                                [20.999999999999993, 0.09962264150943395],
                                [29.39999999999999, 0.0820125786163522],
                                [41, 0.06490566037735851],
                                [55, 0.052327044025157254],
                                [70.19999999999999, 0.04301886792452836],
                                [87.4, 0.03622641509433966],
                                [107.40000000000002, 0.030188679245283068],
                                [133.40000000000003, 0.02515723270440262],
                                [156.60000000000002, 0.02163522012578628],
                                [175.40000000000003, 0.01962264150943399],
                                [200.20000000000005, 0.017358490566037776]])
                from scipy.interpolate import interp1d
                qmin_interp = interp1d(dat[:, 0], dat[:, 1])
                qmin_msort = np.ones_like(primary_mass_high) * 0.1
                ind_5, = np.where(primary_mass_high > 5.0)
                qmin_msort[ind_5] = qmin_interp(primary_mass_high[ind_5])

                q_msort = np.random.uniform(qmin_msort, 1.0)
                secondary_mass_high = q_msort * primary_mass_high

        secondary_mass = np.ones(primary_mass.size)
        secondary_mass[mlowIdx] = secondary_mass_low
        if msort < 10000.0:
            secondary_mass[mhighIdx] = secondary_mass_high

        return secondary_mass

    def binary_select(self, primary_mass, binfrac_model=0.5, **kwargs):
        """Select which primary masses will have a companion using
        either a binary fraction specified by a float or a
        primary-mass dependent binary fraction following
        `van Haaften et al.(2009) <http://adsabs.harvard.edu/abs/2013A%26A...552A..69V>`_ in appdx

        Parameters
        ----------
            primary_mass : array
                Mass that determines the binary fraction

            binfrac_model : str or float
                vanHaaften - primary mass dependent and ONLY VALID up to 100 Msun
                float - fraction of binaries; 0.5 means 2 in 3 stars are a binary pair while 1
                means every star is in a binary pair

        Returns
        -------
            stars_in_binary : array
                primary masses that will have a binary companion

            stars_in_single : array
                primary masses that will be single stars

            binary_fraction : array
                system-specific probability of being in a binary

            binaryIdx : array
                Idx of stars in binary
        """

        flag_msort = kwargs.pop("flag_msort", "no_msort")

        if flag_msort == "binfrac_high_mass":
            msort = kwargs.pop("msort")
            binfrac_model_msort = kwargs.pop("binfrac_model_msort")
        else:
            msort, binfrac_model_msort = 10000., 0.

        sorting = msort * np.ones(primary_mass.size)
        (mhighIdx,) = np.where(primary_mass >= sorting)
        (mlowIdx,) = np.where(primary_mass < sorting)
        primary_mass_high = primary_mass[mhighIdx]
        primary_mass_low = primary_mass[mlowIdx]

        if type(binfrac_model) == str:
            if binfrac_model == "vanHaaften":
                binary_fraction_low = 1 / 2.0 + 1 / \
                    4.0 * np.log10(primary_mass_low)
                binary_choose_low = np.random.uniform(
                    0, 1.0, primary_mass_low.size)

                (singleIdx_low,) = np.where(
                    binary_fraction_low < binary_choose_low)
                (binaryIdx_low,) = np.where(
                    binary_fraction_low >= binary_choose_low)
            else:
                raise ValueError(
                    "You have supplied a non-supported binary fraction model. Please choose vanHaaften or a float"
                )
        elif type(binfrac_model) == float:
            if (binfrac_model <= 1.0) & (binfrac_model >= 0.0):
                binary_fraction_low = binfrac_model * \
                    np.ones(primary_mass_low.size)
                binary_choose_low = np.random.uniform(
                    0, 1.0, primary_mass_low.size)

                (singleIdx_low,) = np.where(
                    binary_choose_low > binary_fraction_low)
                (binaryIdx_low,) = np.where(
                    binary_choose_low <= binary_fraction_low)
            else:
                raise ValueError(
                    "You have supplied a fraction outside of 0-1. Please choose a fraction between 0 and 1."
                )
        else:
            raise ValueError(
                "You have not supplied a model or a fraction. Please choose either vanHaaften or a float"
            )

        if msort < 10000. and type(binfrac_model_msort) == str:
            if binfrac_model_msort == "vanHaaften":
                binary_fraction_high = 1 / 2.0 + 1 / \
                    4.0 * np.log10(primary_mass_high)
                binary_choose_high = np.random.uniform(
                    0, 1.0, primary_mass_high.size)

                (singleIdx_high,) = np.where(
                    binary_fraction_high < binary_choose_high)
                (binaryIdx_high,) = np.where(
                    binary_fraction_high >= binary_choose_high)
            else:
                raise ValueError(
                    "You have supplied a non-supported binary fraction model. Please choose vanHaaften or a float"
                )
        elif msort < 10000. and type(binfrac_model_msort) == float:
            if (binfrac_model_msort <= 1.0) & (binfrac_model_msort >= 0.0):
                binary_fraction_high = binfrac_model_msort * \
                    np.ones(primary_mass_high.size)
                binary_choose_high = np.random.uniform(
                    0, 1.0, primary_mass_high.size)

                (singleIdx_high,) = np.where(
                    binary_choose_high > binary_fraction_high)
                (binaryIdx_high,) = np.where(
                    binary_choose_high <= binary_fraction_high)
            else:
                raise ValueError(
                    "You have supplied a fraction outside of 0-1. Please choose a fraction between 0 and 1."
                )
        elif msort < 10000.:
            raise ValueError(
                "You have not supplied a model or a fraction. Please choose either vanHaaften or a float"
            )

        if msort < 10000.:
            stars_in_binary = np.append(
                primary_mass_high[binaryIdx_high], primary_mass_low[binaryIdx_low])
            stars_in_single = np.append(
                primary_mass_high[singleIdx_high], primary_mass_low[singleIdx_low])
            binary_fraction = np.append(
                binary_fraction_high[binaryIdx_high], binary_fraction_low[binaryIdx_low])
            binaryIdx = np.append(
                mhighIdx[binaryIdx_high], mlowIdx[binaryIdx_low])
        else:
            stars_in_binary = primary_mass_low[binaryIdx_low]
            stars_in_single = primary_mass_low[singleIdx_low]
            binary_fraction = binary_fraction_low[binaryIdx_low]
            binaryIdx = mlowIdx[binaryIdx_low]

        return (
            stars_in_binary,
            stars_in_single,
            binary_fraction,
            binaryIdx,
        )

    def sample_porb(self, mass1, mass2, rad1, rad2, porb_model="sana12", porb_max=None, size=None):
        """Sample the orbital period according to the user-specified model

        Parameters
        ----------
        mass1 : array
            primary masses
        mass2 : array
            secondary masses
        rad1 : array
            radii of the primaries. 
        rad2 : array
            radii of the secondaries 
        model : string
            selects which model to sample orbital periods, choices include:
            log_uniform : semi-major axis flat in log space from RRLO < 0.5 up to 1e5 Rsun according to
            `Abt (1983) <http://adsabs.harvard.edu/abs/1983ARA%26A..21..343A>`_
            and consistent with Dominik+2012,2013
            and then converted to orbital period in days using Kepler III
            sana12 : power law orbital period between 0.15 < log(P/day) < 5.5 following
            `Sana+2012 <https://ui.adsabs.harvard.edu/abs/2012Sci...337..444S/abstract>_`
            renzo19 : power law orbital period for m1 > 15Msun binaries from
            `Sana+2012 <https://ui.adsabs.harvard.edu/abs/2012Sci...337..444S/abstract>_`
            following the implementation of
            `Renzo+2019 <https://ui.adsabs.harvard.edu/abs/2019A%26A...624A..66R/abstract>_`
            and flat in log otherwise

        Returns
        -------
        porb : array
            orbital period with array size equalling array size
            of mass1 and mass2 in units of days
        aRL_over_a: array
            ratio of radius where RL overflow starts to the sampled seperation
            used to truncate the eccentricitiy distribution
        """

        # First we need to compute where RL overflow starts.  We truncate the lower-bound
        # of the period distribution there
        q = mass2 / mass1
        RL_fac = (0.49 * q ** (2.0 / 3.0)) / (
            0.6 * q ** (2.0 / 3.0) + np.log(1 + q ** 1.0 / 3.0)
        )

        q2 = mass1 / mass2
        RL_fac2 = (0.49 * q2 ** (2.0 / 3.0)) / (
            0.6 * q2 ** (2.0 / 3.0) + np.log(1 + q2 ** 1.0 / 3.0)
        )

        # include the factor for the eccentricity
        RL_max = 2 * rad1 / RL_fac
        (ind_switch,) = np.where(RL_max < 2 * rad2 / RL_fac2)
        if len(ind_switch) >= 1:
            RL_max[ind_switch] = 2 * rad2[ind_switch] / RL_fac2[ind_switch]

        # Can either sample the porb first and truncate the eccentricities at RL overflow
        # or sample the eccentricities first and truncate a(1-e) at RL overflow
        #
        # If we haven't sampled the eccentricities, then the minimum semi-major axis is at
        # RL overflow
        #
        # If we have, then the minimum pericenter is set to RL overflow
        a_min = RL_max 

        if porb_model == "log_uniform":
            if porb_max is None:
                a_0 = np.random.uniform(np.log(a_min), np.log(1e5), size)
            else:
                # If in CMC, only sample binaries as wide as the local hard/soft boundary
                a_max = utils.a_from_p(porb_max,mass1,mass2) 
                a_max[a_max < a_min] = a_min[a_max < a_min]
                a_0 = np.random.uniform(np.log(a_min), np.log(a_max), size)

            # convert out of log space
            a_0 = np.exp(a_0)
            aRL_over_a = a_min/a_0

            # convert to au
            rsun_au = 0.00465047
            a_0 = a_0 * rsun_au

            # convert to orbital period in years
            yr_day = 365.24
            porb_yr = ((a_0 ** 3.0) / (mass1 + mass2)) ** 0.5
            porb = porb_yr * yr_day
        elif porb_model == "sana12":
            # Same here: if using CMC, set the maximum porb to the smaller of either the
            # hard/soft boundary or 5.5 (from Sana paper)
            if porb_max is None:
                log10_porb_max = 5.5
            else:
                log10_porb_max = np.minimum(5.5,np.log10(porb_max))

            porb = 10 ** utils.rndm(a=0.15, b=log10_porb_max, g=-0.55, size=size)
            aRL_over_a = a_min / utils.a_from_p(porb,mass1,mass2) 

        elif porb_model == "renzo19":
            # Same here: if using CMC, set the maximum porb to the smaller of either the
            # hard/soft boundary or 5.5 (from Sana paper)
            if porb_max is None:
                log10_porb_max = 5.5
            else:
                log10_porb_max = np.minimum(5.5,np.log10(porb_max))

            porb = 10 ** (np.random.uniform(0.15, log10_porb_max, size))
            (ind_massive,) = np.where(mass1 > 15)
            porb[ind_massive] = 10 ** utils.rndm(
                a=0.15, b=log10_porb_max, g=-0.55, size=len(ind_massive)
            )
            aRL_over_a = a_min / utils.a_from_p(porb,mass1,mass2) 
        else:
            raise ValueError(
                "You have supplied a non-supported model; Please choose either log_flat, sana12, or renzo19"
            )
        return porb, aRL_over_a 

    def sample_ecc(self, aRL_over_a, ecc_model="sana12", size=None):
        """Sample the eccentricity according to a user specified model

        Parameters
        ----------
        ecc_model : string
            'thermal' samples from a  thermal eccentricity distribution following
            `Heggie (1975) <http://adsabs.harvard.edu/abs/1975MNRAS.173..729H>`_
            'uniform' samples from a uniform eccentricity distribution
            'sana12' samples from the eccentricity distribution from
            `Sana+2012 <https://ui.adsabs.harvard.edu/abs/2012Sci...337..444S/abstract>_`
            'circular' assumes zero eccentricity for all systems
            DEFAULT = 'sana12'

        aRL_over_a : ratio of the minimum seperation (where RL overflow starts)
            to the sampled semi-major axis.  Use this to truncate the eccentricitiy

        size : int, optional
            number of eccentricities to sample
            this is set in cosmic-pop call as Nstep

        Returns
        -------
        ecc : array
            array of sampled eccentricities with size=size
        """

        # if we sampled the periods first, we need to truncate the eccentricities
        # to avoid RL overflow/collision at pericenter
        e_max = 1.0 - aRL_over_a

        if ecc_model == "thermal":
            a_0 = np.random.uniform(0.0, e_max**2, size)
            ecc = a_0 ** 0.5
            return ecc

        elif ecc_model == "uniform":
            ecc = np.random.uniform(0.0, e_max, size)
            return ecc

        elif ecc_model == "sana12":
            ecc = utils.rndm(a=0.001, b=0.9, g=-0.45, size=size)
            return ecc

        elif ecc_model == "circular":
            ecc = np.zeros(size)
            return ecc

        else:
            raise ValueError("You have specified an unsupported model. Please choose from thermal, "
                             "uniform, sana12, or circular")

    def sample_SFH(self, SF_start=13700.0, SF_duration=0.0, met=0.02, size=None):
        """Sample an evolution time for each binary based on a user-specified
        time at the start of star formation and the duration of star formation.
        The default is a burst of star formation 13,700 Myr in the past.

        Parameters
        ----------
        SF_start : float
            Time in the past when star formation initiates in Myr
        SF_duration : float
            Duration of constant star formation beginning from SF_Start in Myr
        met : float
            metallicity of the population [Z_sun = 0.02]
            Default: 0.02
        size : int, optional
            number of evolution times to sample
            NOTE: this is set in cosmic-pop call as Nstep

        Returns
        -------
        tphys : array
            array of evolution times of size=size
        metallicity : array
            array of metallicities
        """

        if (SF_start > 0.0) & (SF_duration >= 0.0):
            tphys = np.random.uniform(SF_start - SF_duration, SF_start, size)
            metallicity = np.ones(size) * met
            return tphys, metallicity
        else:
            raise ValueError(
                'SF_start and SF_duration must be positive and SF_start must be greater than 0.0')

    def set_kstar(self, mass):
        """Initialize stellar types according to BSE classification
        kstar=1 if M>=0.7 Msun; kstar=0 if M<0.7 Msun

        Parameters
        ----------
        mass : array
            array of masses

        Returns
        -------
        kstar : array
            array of initial stellar types
        """

        kstar = np.zeros(mass.size)
        low_cutoff = 0.7
        lowIdx = np.where(mass < low_cutoff)[0]
        hiIdx = np.where(mass >= low_cutoff)[0]

        kstar[lowIdx] = 0
        kstar[hiIdx] = 1

        return kstar

    def set_reff(self, mass, metallicity, zsun=0.02):
        """
        Better way to set the radii from BSE, by calling it directly

        takes masses and metallicities, and returns the radii

        Note that the BSE function is hard-coded to go through arrays
        of length 10^5.  If your masses are more than that, you'll
        need to divide it into chunks
        """

        max_array_size = 100000
        total_length = len(mass)
        radii = np.zeros(total_length)

        _evolvebin.metvars.zsun = zsun

        idx = 0
        while total_length > max_array_size:
            ## cycle through the masses max_array_size number at a time
            temp_mass = mass[idx*max_array_size:(idx+1)*max_array_size]

            temp_radii = _evolvebin.compute_r(temp_mass,metallicity,max_array_size)

            ## put these in the radii array
            radii[idx*max_array_size:(idx+1)*max_array_size] = temp_radii

            total_length -= max_array_size
            idx += 1

        length_remaining = total_length

        ## if smaller than 10^5, need to pad out the array
        temp_mass = np.zeros(max_array_size)
        temp_mass[:length_remaining] = mass[-length_remaining:]

        temp_radii = _evolvebin.compute_r(temp_mass,metallicity,length_remaining)

        #finish up the array
        radii[-length_remaining:] = temp_radii[:length_remaining]

        return radii
