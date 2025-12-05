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
import warnings
import pandas as pd

from cosmic import utils

from .sampler import register_sampler
from .. import InitialBinaryTable


__author__ = "Katelyn Breivik <katie.breivik@gmail.com>"
__credits__ = ("Scott Coughlin <scott.coughlin@ligo.org>, Michael Zevin <michael.j.zevin@gmail.com>, "
               "Tom Wagg <tomjwagg@gmail.com>")
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
    size=None,
    total_mass=np.inf,
    sampling_target="size",
    trim_extra_samples=False,
    q_power_law=0,
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

    porb_model : `str` or `dict`
        Model to sample orbital period; choices include: log_uniform, sana12, raghavan10, moe19
        or a custom power law distribution defined with a dictionary with keys "min", "max", and "slope"
        (e.g. {"min": 0.15, "max": 0.55, "slope": -0.55}) would reproduce the Sana+2012 distribution

    qmin : `float`
        kwarg which sets the minimum mass ratio for sampling the secondary
        where the mass ratio distribution is flat in q
        if q > 0, qmin sets the minimum mass ratio
        q = -1, this limits the minimum mass ratio to be set such that
        the pre-MS lifetime of the secondary is not longer than the full
        lifetime of the primary if it were to evolve as a single star

    m_max : `float`
        kwarg which sets the maximum primary and secondary mass for sampling
        NOTE: this value changes the range of the IMF and should *not* be used
            as a means of selecting certain kstar types!

    m1_min : `float`
        kwarg which sets the minimum primary mass for sampling
        NOTE: this value changes the range of the IMF and should *not* be used
            as a means of selecting certain kstar types!

    m2_min : `float`
        kwarg which sets the minimum secondary mass for sampling
        the secondary as uniform in mass_2 between m2_min and mass_1

    msort : `float`
        Stars with M>msort can have different pairing and sampling of companions

    qmin_msort : `float`
        Same as qmin for M>msort

    m2_min_msort : `float`
        Same as m2_min for M>msort

    SF_start : `float`
        Time in the past when star formation initiates in Myr

    SF_duration : `float`
        Duration of constant star formation beginning from SF_Start in Myr

    binfrac_model : `str or float`
        Model for binary fraction; choices include: vanHaaften, offner22, or a fraction where 1.0 is 100% binaries

    binfrac_model_msort : `str or float`
        Same as binfrac_model for M>msort

    met : `float`
        Sets the metallicity of the binary population where solar metallicity is zsun

    size : `int`
        Size of the population to sample

    total_mass : `float`
        Total mass to use as a target for sampling

    sampling_target : `str`
        Which type of target to use for sampling (either "size" or "total_mass"), by default "size".
        Note that `total_mass` must not be None when `sampling_target=="total_mass"`.

    trim_extra_samples : `str`
        Whether to trim the sampled population so that the total mass sampled is as close as possible to
        `total_mass`. Ignored when `sampling_target==size`.
        Note that given the discrete mass of stars, this could mean your sample is off by 300
        solar masses in the worst case scenario (of a 150+150 binary being sampled). In reality the majority
        of cases track the target total mass to within a solar mass.        

    zsun : `float`
        optional kwarg for setting effective radii, default is 0.02

    q_power_law : `float`
        Exponent for the mass ratio distribution power law, default is 0 (flat in q). Note that
        q_power_law cannot be exactly -1, as this would result in a divergent distribution.


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
    if sampling_target == "total_mass" and (total_mass is None or total_mass == np.inf):
        raise ValueError("If `sampling_target == 'total mass'` then `total_mass` must be supplied")
    if size is None and (total_mass is None or total_mass == np.inf):
        raise ValueError("Either a sample `size` or `total_mass` must be supplied")
    elif size is None:
        size = int(total_mass)

    if binfrac_model == 0.0 and sampling_target == "size":
        raise ValueError(("If `binfrac_model == 0.0` then `sampling_target` must be 'total_mass'. Otherwise "
                          "you are targetting a population of `size` binaries but will never select any."))

    final_kstar1 = [final_kstar1] if isinstance(final_kstar1, (int, float)) else final_kstar1
    final_kstar2 = [final_kstar2] if isinstance(final_kstar2, (int, float)) else final_kstar2
    primary_min, primary_max, secondary_min, secondary_max = utils.mass_min_max_select(
        final_kstar1, final_kstar2, **kwargs)
    initconditions = Sample()

    # set up multiplier if the mass sampling is inefficient
    multiplier = 1

    # track samples to actually return (after masks)
    mass1_singles = []
    mass1_binary = []
    mass2_binary = []
    binfrac = []

    # track the total mass of singles and binaries sampled
    m_sampled_singles = 0.0
    m_sampled_binaries = 0.0

    # track the total number of stars sampled
    n_singles = 0
    n_binaries = 0

    # if porb_model = `moe19`, the binary fraction is fixed based on the metallicity 
    if porb_model == "moe19":
        binfrac_model = utils.get_met_dep_binfrac(met)
        warnings.warn('your supplied binfrac_model has been overwritten to {} match Moe+2019'.format(binfrac_model))

    # define a function that evaluates whether you've reached your sampling target
    target = lambda mass1_binary, size, m_sampled_singles, m_sampled_binaries, total_mass:\
        len(mass1_binary) < size if sampling_target == "size" else m_sampled_singles + m_sampled_binaries < total_mass

    # sample until you've reached your target
    while target(mass1_binary, size, m_sampled_singles, m_sampled_binaries, total_mass):
        # sample primary masses
        mass1, _ = initconditions.sample_primary(primary_model, size=int(size * multiplier), **kwargs)

        # split them into binaries or single stars
        (mass1_binaries, mass_single, binfrac_binaries, binary_index,
        ) = initconditions.binary_select(mass1, binfrac_model=binfrac_model, **kwargs)

        # sample secondary masses for the single stars
        mass2_binaries = initconditions.sample_secondary(mass1_binaries, q_power_law=q_power_law, **kwargs)

        # check if this batch of samples will take us over our sampling target
        if not target(mass1_binary, size,
                      m_sampled_singles + np.sum(mass_single),
                      m_sampled_binaries + np.sum(mass1_binaries) + np.sum(mass2_binaries),
                      total_mass) and trim_extra_samples and sampling_target == "total_mass":
            # get the cumulative total mass of the samples
            total_mass_list = np.copy(mass1)
            total_mass_list[binary_index] += mass2_binaries
            sampled_so_far = m_sampled_singles + m_sampled_binaries
            cumulative_total_mass = sampled_so_far + np.cumsum(total_mass_list)

            # find the boundary for reaching the right total mass
            threshold_index = np.where(cumulative_total_mass > total_mass)[0][0]

            keep_offset = abs(cumulative_total_mass[threshold_index] - total_mass)
            drop_offset = abs(cumulative_total_mass[threshold_index - 1] - total_mass)
            lim = threshold_index - 1 if (keep_offset > drop_offset) else threshold_index
            
            # work out how many singles vs. binaries to delete
            one_if_binary = np.zeros(len(mass1))
            one_if_binary[binary_index] = 1
            sb_delete = one_if_binary[lim + 1:]
            n_single_delete = (sb_delete == 0).sum()
            n_binary_delete = (sb_delete == 1).sum()

            # delete em!
            if n_single_delete > 0:
                mass_single = mass_single[:-n_single_delete]
            if n_binary_delete > 0:
                mass1_binaries = mass1_binaries[:-n_binary_delete]
                mass2_binaries = mass2_binaries[:-n_binary_delete]
                binfrac_binaries = binfrac_binaries[:-n_binary_delete]

            # ensure we don't loop again after this
            target = lambda mass1_binary, size, m_sampled_singles, m_sampled_binaries, total_mass: False

        # track the mass sampled
        m_sampled_singles += sum(mass_single)
        m_sampled_binaries += sum(mass1_binaries)
        m_sampled_binaries += sum(mass2_binaries)

        # track the total number sampled
        n_singles += len(mass_single)
        n_binaries += len(mass1_binaries)

        # select out the primaries and secondaries that will produce the final kstars
        ind_select = (  (mass1_binaries > primary_min)
                      & (mass1_binaries < primary_max)
                      & (mass2_binaries > secondary_min)
                      & (mass2_binaries < secondary_max))
        mass1_binary.extend(mass1_binaries[ind_select])
        mass2_binary.extend(mass2_binaries[ind_select])
        binfrac.extend(binfrac_binaries[ind_select])

        # select out the single stars that will produce the final kstar
        mass1_singles.extend(mass_single[(mass_single > primary_min) & (mass_single < primary_max)])

        # check to see if we should increase the multiplier factor to sample the population more quickly
        if target(mass1_binary, size / 100, m_sampled_singles, m_sampled_binaries, total_mass / 100):
            # well this sampling rate is clearly not working time to increase
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
    # if the porb_model is moe19, the metallicity needs to be supplied
    if porb_model == "moe19":
        porb,aRL_over_a = initconditions.sample_porb(
            mass1_binary, mass2_binary, rad1, rad2, porb_model, met=met, size=mass1_binary.size
        )
    else:
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
        tphysf_singles, metallicity_singles = initconditions.sample_SFH(
            SF_start=SF_start, SF_duration=SF_duration, met=met, size=mass1_singles.size
        )
        metallicity_singles[metallicity_singles < 1e-4] = 1e-4
        metallicity_singles[metallicity_singles > 0.03] = 0.03
        kstar1_singles = initconditions.set_kstar(mass1_singles)
        singles_table = InitialBinaryTable.InitialBinaries(
            mass1_singles,                          # mass1
            np.ones_like(mass1_singles) * 0,        # mass2 (all massless remnants)
            np.ones_like(mass1_singles) * -1,       # porb (single not binary)
            np.ones_like(mass1_singles) * -1,       # ecc (single not binary)
            tphysf_singles,                         # tphysf
            kstar1_singles,                         # kstar1
            np.ones_like(mass1_singles) * 15,       # kstar2 (all massless remnants)
            metallicity_singles,                    # metallicity
        )
        binary_table = pd.concat([binary_table, singles_table], ignore_index=True)
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
        m_sampled_singles,
        m_sampled_binaries,
        n_singles,
        n_binaries
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

            Optional kwargs are defined in `get_independent_sampler`

            Returns
            -------
            a_0 : array
            Sampled primary masses
            np.sum(a_0) : float
            Total amount of mass sampled
            """

        # Read in m1_min and m_max kwargs, if provided
        m1_min = kwargs["m1_min"] if "m1_min" in kwargs.keys() else 0.08
        m_max = kwargs["m_max"] if "m_max" in kwargs.keys() else 150.0

        # Make sure m1_min value is below 0.5, since otherwise it will not work for Kroupa IMF
        if m1_min > 0.5:
            raise ValueError("m1_min must be greater than 0.5 Msun")

        if primary_model == 'kroupa93':
            alphas, mcuts = [-1.3,-2.2,-2.7], [m1_min,0.5,1.0,m_max]
        # Since COSMIC/BSE can't handle < 0.08Msun, we by default truncate at 0.08 Msun instead of 0.01
        elif primary_model == 'kroupa01':
            alphas, mcuts = [-1.3,-2.3], [m1_min,0.5,m_max]
        elif primary_model == 'salpeter55':
            alphas, mcuts = [-2.35], [m1_min,m_max]
        elif primary_model == 'custom':
            if 'alphas' in kwargs and 'mcuts' in kwargs:
                alphas = kwargs.pop("alphas", [-1.3,-2.3,-2.3])
                mcuts = kwargs.pop("mcuts", [m1_min,0.5,1.0,m_max])
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
    def sample_secondary(self, primary_mass, q_power_law=0, **kwargs):
        """Sample a secondary mass using draws from a uniform mass ratio distribution motivated by
        `Mazeh et al. (1992) <http://adsabs.harvard.edu/abs/1992ApJ...401..265M>`_
        and `Goldberg & Mazeh (1994) <http://adsabs.harvard.edu/abs/1994ApJ...429..362G>`_

        NOTE: the lower lim is set by either qmin or m2_min which are passed as kwargs

        Parameters
        ----------
        primary_mass : `array`
            sets the maximum secondary mass (for a maximum mass ratio of 1)

        Optional kwargs are defined in `get_independent_sampler`

        Returns
        -------
        secondary_mass : array
            sampled secondary masses with array size matching size of
            primary_mass
        """
        qmin = kwargs["qmin"] if "qmin" in kwargs.keys() else 0.0
        m1_min = kwargs["m1_min"] if "m1_min" in kwargs.keys() else 0.08
        m2_min = kwargs["m2_min"] if "m2_min" in kwargs.keys() else None
        if (m2_min is None) & (qmin is None):
            warnings.warn("It is highly recommended that you specify either qmin or m2_min!")
        if (m2_min is not None) and (m2_min > m1_min):
            raise ValueError("The m2_min you specified is above the minimum"
                             " primary mass of the IMF, either lower m2_min or"
                             " raise the lower value of your sampled primaries")
        
        if (m2_min is not None) & (qmin != 0):
            raise ValueError("You cannot specify both m2_min and qmin, please choose one or the other")  

        # --- `msort` kwarg can be set to have different qmin above `msort`
        msort = kwargs["msort"] if "msort" in kwargs.keys() else None
        qmin_msort = kwargs["qmin_msort"] if "qmin_msort" in kwargs.keys() else None
        m2_min_msort = kwargs["m2_min_msort"] if "m2_min_msort" in kwargs.keys() else None
        if (msort is None) and (qmin_msort is not None):
            raise ValueError("If qmin_msort is specified, you must also supply a value for msort")
        if (msort is None) and (m2_min_msort is not None):
            raise ValueError("If m2_min_msort is specified, you must also supply a value for msort")
        if (m2_min_msort is not None) and (m2_min_msort > msort):
            raise ValueError("The m2_min_msort you specified is above the minimum"
                             " primary mass of the high-mass binaries msort")

        if (msort is not None) and (qmin_msort is not None):
            (highmassIdx,) = np.where(primary_mass >= msort)
            (lowmassIdx,) = np.where(primary_mass < msort)
        else:
            (highmassIdx,) = np.where(primary_mass < 0)
            (lowmassIdx,) = np.where(primary_mass >= 0)   # all idxs


        qmin_vals = -1 * np.ones_like(primary_mass)
        # --- qmin for low-mass systems (all systems if msort is not specified)
        if (qmin > 0.0):
            qmin_vals[lowmassIdx] = qmin * np.ones_like(primary_mass[lowmassIdx])
        elif (qmin < 0.0):
            # mass-dependent qmin, assume qmin=0.1 for m_primary<5
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
            qmin_vals[lowmassIdx] = np.ones_like(primary_mass[lowmassIdx]) * 0.1
            ind_5, = np.where(primary_mass[lowmassIdx] > 5.0)
            qmin_vals[lowmassIdx][ind_5] = qmin_interp(primary_mass[lowmassIdx][ind_5])
        else:
            qmin_vals[lowmassIdx] = np.zeros_like(primary_mass[lowmassIdx])
        # --- qmin for high-mass systems, if msort and qmin_msort are specified
        if (msort is not None) and (qmin_msort is not None):
            if (qmin_msort > 0.0):
                qmin_vals[highmassIdx] = qmin_msort * np.ones_like(primary_mass[highmassIdx])
            elif (qmin_msort < 0.0):
                # mass-dependent qmin, assume qmin=0.1 for m_primary<5
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
                qmin_vals[highmassIdx] = np.ones_like(primary_mass[highmassIdx]) * 0.1
                ind_5, = np.where(primary_mass[highmassIdx] > 5.0)
                qmin_vals[highmassIdx][ind_5] = qmin_interp(primary_mass[highmassIdx][ind_5])
            else:
                qmin_vals[highmassIdx] = np.zeros_like(primary_mass[highmassIdx])

        # --- apply m2_min and m2_min_msort, if specified
        if m2_min is not None:
            qmin_vals[lowmassIdx] = np.maximum(qmin_vals[lowmassIdx], m2_min/primary_mass[lowmassIdx])
        if m2_min_msort is not None:
            qmin_vals[highmassIdx] = np.maximum(qmin_vals[highmassIdx], m2_min_msort/primary_mass[highmassIdx])

        # --- now, randomly sample mass ratios and get secondary masses
        secondary_mass = utils.rndm(qmin_vals, 1, q_power_law, size=len(primary_mass)) * primary_mass
        return secondary_mass

    def binary_select(self, primary_mass, binfrac_model=0.5, **kwargs):
        """Select which primary masses will have a companion using
        either a binary fraction specified by a float or a
        primary-mass dependent binary fraction following
        `van Haaften et al.(2009) <http://adsabs.harvard.edu/abs/2013A%26A...552A..69V>`_ in appdx
        or `Offner et al.(2022) <https://arxiv.org/abs/2203.10066>`_ in fig 1

        Parameters
        ----------
            primary_mass : array
                Mass that determines the binary fraction

            binfrac_model : str or float
                vanHaaften - primary mass dependent and ONLY VALID up to 100 Msun
                
                offner22 - primary mass dependent
                
                float - fraction of binaries; 0.5 means 2 in 3 stars are a binary pair while 1
                means every star is in a binary pair

        Optional kwargs are defined in `get_independent_sampler`

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

        # --- `msort` kwarg can be set to have different binary fraction above `msort`
        msort = kwargs["msort"] if "msort" in kwargs.keys() else None
        binfrac_model_msort = kwargs["binfrac_model_msort"] if "binfrac_model_msort" in kwargs.keys() else None
        if (msort is None) and (binfrac_model_msort is not None):
            raise ValueError("If binfrac_model_msort is specified, you must also supply a value for msort")
        if (msort is not None) and (binfrac_model_msort is not None):
            (highmassIdx,) = np.where(primary_mass >= msort)
            (lowmassIdx,) = np.where(primary_mass < msort)
        else:
            (highmassIdx,) = np.where(primary_mass < 0)
            (lowmassIdx,) = np.where(primary_mass >= 0)   # all idxs


        # --- read in binfrac models
        if type(binfrac_model) == str:
            if binfrac_model == "vanHaaften":
                binary_fraction_low = 1 / 2.0 + 1 / \
                    4.0 * np.log10(primary_mass[lowmassIdx])
                binary_choose_low = np.random.uniform(
                    0, 1.0, primary_mass[lowmassIdx].size)

                (singleIdx_low,) = np.where(
                    binary_fraction_low < binary_choose_low)
                (binaryIdx_low,) = np.where(
                    binary_fraction_low >= binary_choose_low)
            elif binfrac_model == "offner22":
                from scipy.interpolate import BSpline
                t = [0.0331963853, 0.0331963853, 0.0331963853, 0.0331963853, 0.106066017,
                     0.212132034, 0.424264069, 0.866025404, 1.03077641, 1.11803399,
                     1.95959179, 3.87298335, 6.32455532, 11.6619038, 29.1547595,
                     40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 150, 150, 150, 150]
                c = [0.08, 0.15812003, 0.20314101, 0.23842953, 0.33154153, 0.39131739,
                     0.46020725, 0.59009569, 0.75306454, 0.81652502, 0.93518422, 0.92030594,
                     0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96]
                k = 3
                def offner_curve(x):
                    a = -0.16465041
                    b = -0.11616329
                    return np.piecewise(x, [x < 6.4, x >= 6.4], [BSpline(t,c,k), lambda x : a * np.exp(b * x) + 0.97])
                binary_fraction_low = offner_curve(primary_mass[lowmassIdx])
                binary_choose_low = np.random.uniform(
                    0, 1.0, primary_mass[lowmassIdx].size)

                (singleIdx_low,) = np.where(
                    binary_fraction_low < binary_choose_low)
                (binaryIdx_low,) = np.where(
                    binary_fraction_low >= binary_choose_low)
            else:
                raise ValueError(
                    "You have supplied a non-supported binary fraction model. Please choose vanHaaften, offner22, or a float"
                )
        elif type(binfrac_model) == float:
            if (binfrac_model <= 1.0) & (binfrac_model >= 0.0):
                binary_fraction_low = binfrac_model * \
                    np.ones(primary_mass[lowmassIdx].size)
                binary_choose_low = np.random.uniform(
                    0, 1.0, primary_mass[lowmassIdx].size)

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
                "You have not supplied a model or a fraction. Please choose either vanHaaften, offner22, or a float"
            )

        # --- if using a different binary fraction for high-mass systems
        if (binfrac_model_msort is not None) and (type(binfrac_model_msort) == str):
            if binfrac_model_msort == "vanHaaften":
                binary_fraction_high = 1 / 2.0 + 1 / \
                    4.0 * np.log10(primary_mass[highmassIdx])
                binary_choose_high = np.random.uniform(
                    0, 1.0, primary_mass[highmassIdx].size)

                (singleIdx_high,) = np.where(
                    binary_fraction_high < binary_choose_high)
                (binaryIdx_high,) = np.where(
                    binary_fraction_high >= binary_choose_high)
            elif binfrac_model_msort == "offner22":
                from scipy.interpolate import BSpline
                t = [0.0331963853, 0.0331963853, 0.0331963853, 0.0331963853, 0.106066017,
                     0.212132034, 0.424264069, 0.866025404, 1.03077641, 1.11803399,
                     1.95959179, 3.87298335, 6.32455532, 11.6619038, 29.1547595,
                     40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 150, 150, 150, 150]
                c = [0.08, 0.15812003, 0.20314101, 0.23842953, 0.33154153, 0.39131739,
                     0.46020725, 0.59009569, 0.75306454, 0.81652502, 0.93518422, 0.92030594,
                     0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96]
                k = 3
                def offner_curve(x):
                    a = -0.16465041
                    b = -0.11616329
                    return np.piecewise(x, [x < 6.4, x >= 6.4], [BSpline(t,c,k), lambda x : a * np.exp(b * x) + 0.97])
                binary_fraction_high = offner_curve(primary_mass[highmassIdx])
                binary_choose_high = np.random.uniform(
                    0, 1.0, primary_mass[highmassIdx].size)

                (singleIdx_high,) = np.where(
                    binary_fraction_high < binary_choose_high)
                (binaryIdx_high,) = np.where(
                    binary_fraction_high >= binary_choose_high)
            else:
                raise ValueError(
                    "You have supplied a non-supported binary fraction model. Please choose vanHaaften, offner22, or a float"
                )
        elif (binfrac_model_msort is not None) and (type(binfrac_model_msort) == float):
            if (binfrac_model_msort <= 1.0) & (binfrac_model_msort >= 0.0):
                binary_fraction_high = binfrac_model_msort * \
                    np.ones(primary_mass[highmassIdx].size)
                binary_choose_high = np.random.uniform(
                    0, 1.0, primary_mass[highmassIdx].size)

                (singleIdx_high,) = np.where(
                    binary_choose_high > binary_fraction_high)
                (binaryIdx_high,) = np.where(
                    binary_choose_high <= binary_fraction_high)
            else:
                raise ValueError(
                    "You have supplied a fraction outside of 0-1. Please choose a fraction between 0 and 1."
                )
        elif (binfrac_model_msort is not None):
            raise ValueError(
                "You have not supplied a model or a fraction. Please choose either vanHaaften, offner22, or a float"
            )


        # --- get pertinent info
        if (binfrac_model_msort is not None):
            stars_in_binary = np.append(
                primary_mass[highmassIdx][binaryIdx_high], primary_mass[lowmassIdx][binaryIdx_low])
            stars_in_single = np.append(
                primary_mass[highmassIdx][singleIdx_high], primary_mass[lowmassIdx][singleIdx_low])
            binary_fraction = np.append(
                binary_fraction_high[binaryIdx_high], binary_fraction_low[binaryIdx_low])
            binaryIdx = np.append(
                highmassIdx[binaryIdx_high], lowmassIdx[binaryIdx_low])
        else:
            stars_in_binary = primary_mass[lowmassIdx][binaryIdx_low]
            stars_in_single = primary_mass[lowmassIdx][singleIdx_low]
            binary_fraction = binary_fraction_low[binaryIdx_low]
            binaryIdx = lowmassIdx[binaryIdx_low]

        return (
            stars_in_binary,
            stars_in_single,
            binary_fraction,
            binaryIdx,
        )

    def sample_porb(self, mass1, mass2, rad1, rad2, porb_model, porb_max=None, size=None, **kwargs):
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
        porb_model : `str` or `dict`
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
            raghavan10 : log normal orbital periods in days with mean_logP = 4.9 
            and sigma_logP = 2.3 between 0 < log10(P/day) < 9  following 
            `Raghavan+2010 <https://ui.adsabs.harvard.edu/abs/2010ApJS..190....1R/abstract>_`
            moe19 : log normal orbital periods in days with mean_logP = 4.9 
            and sigma_logP = 2.3 between 0 < log10(P/day) < 9  following 
            `Raghavan+2010 <https://ui.adsabs.harvard.edu/abs/2010ApJS..190....1R/abstract>_`
            but with different close binary fractions following 
            `Moe+2019 <https://ui.adsabs.harvard.edu/abs/2019ApJ...875...61M/abstract>_`
            Custom power law distribution defined with a dictionary with keys "min", "max", and "slope"
            (e.g. porb_model={"min": 0.15, "max": 0.55, "slope": -0.55}) would reproduce the
            Sana+2012 distribution.
        met : float
            metallicity of the population

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

            # Use the lower limit from the Sana12 distribution, unless this means the binaries are sampled at RL overflow. If so, 
            # change the lower limit to a_min
                         
            log10_porb_min = np.array([0.15]*len(a_min)) 
            RL_porb = utils.p_from_a(a_min,mass1,mass2)
            log10_RL_porb = np.log10(RL_porb)
            log10_porb_min[log10_porb_min <  log10_RL_porb] = log10_RL_porb[log10_porb_min < log10_RL_porb]
            
            porb = 10 ** utils.rndm(a=log10_porb_min, b=log10_porb_max, g=-0.55, size=size)
            aRL_over_a = a_min / utils.a_from_p(porb,mass1,mass2) 

        elif isinstance(porb_model, dict):
            # use a power law distribution for the orbital periods
            params = {
                "min": 0.15,
                "max": 5.5,
                "slope": -0.55,
            }
            # update the default parameters with the user-supplied ones
            params.update(porb_model)

            # same calculations as sana12 case (sample from a power law distribution but avoid RLOF)
            log10_RL_porb = np.log10(utils.p_from_a(a_min, mass1, mass2))
            params["min"] = np.full(len(a_min), params["min"])
            params["min"][params["min"] < log10_RL_porb] = log10_RL_porb[params["min"] < log10_RL_porb]
            porb = 10**utils.rndm(a=params["min"], b=params["max"], g=params["slope"], size=size)
            aRL_over_a = a_min / utils.a_from_p(porb, mass1, mass2)

        elif porb_model == "renzo19":
            # Same here: if using CMC, set the maximum porb to the smaller of either the
            # hard/soft boundary or 5.5 (from Sana paper)
            if porb_max is None:
                log10_porb_max = 5.5
            
            else:
                log10_porb_max = np.minimum(5.5,np.log10(porb_max))
            
            # Use the lower limit from the Sana12 distribution, unless this means the binaries are sampled at RL overflow. If so, 
            # change the lower limit to a_min
            log10_porb_min = np.array([0.15]*len(a_min)) 
            RL_porb = utils.p_from_a(a_min,mass1,mass2)
            log10_RL_porb = np.log10(RL_porb)
            log10_porb_min[log10_porb_min <  log10_RL_porb] = log10_RL_porb[log10_porb_min < log10_RL_porb]
            
            porb = 10 ** (np.random.uniform(log10_porb_min, log10_porb_max, size))
            (ind_massive,) = np.where(mass1 > 15)
            
            if type(log10_porb_max) != float:
                log10_porb_max = log10_porb_max[ind_massive]
                log10_porb_min = log10_porb_min[ind_massive]              
            
            
            porb[ind_massive] = 10 ** utils.rndm(
                a=log10_porb_min[ind_massive], b=log10_porb_max, g=-0.55, size=len(ind_massive))
            aRL_over_a = a_min / utils.a_from_p(porb,mass1,mass2) 
        
        elif porb_model == "raghavan10":
            import scipy
            # Same here: if using CMC, set the maximum porb to the smaller of either the
            # hard/soft boundary or 5.5 (from Sana paper)
            if porb_max is None:
                log10_porb_max = 9.0
            else:
                log10_porb_max = np.minimum(5.5, np.log10(porb_max))

            lower = 0
            upper = log10_porb_max
            mu = 4.9
            sigma = 2.3

            porb = 10 ** (scipy.stats.truncnorm.rvs(
                (lower-mu)/sigma,(upper-mu)/sigma, loc=mu, scale=sigma, size=size
                ))

            aRL_over_a = a_min / utils.a_from_p(porb,mass1,mass2)

        elif porb_model == "moe19":
            from scipy.interpolate import interp1d
            from scipy.stats import norm
            from scipy.integrate import trapezoid

            try:
                met = kwargs.pop('met')
            except:
                raise ValueError(
                    "You have chosen moe19 for the orbital period distribution which is a metallicity-dependent distribution. "
                    "Please specify a metallicity for the population."
                    )
            def get_logP_dist(nsamp, norm_wide, norm_close, mu=4.4, sigma=2.1):
                logP_lo_lim=0
                logP_hi_lim=9
                close_logP=4.0
                wide_logP=6.0
                neval = 500
                prob_wide = norm.pdf(np.linspace(wide_logP, logP_hi_lim, neval), loc=mu, scale=sigma)*norm_wide
                prob_close = norm.pdf(np.linspace(logP_lo_lim, close_logP, neval), loc=mu, scale=sigma)*norm_close
                slope = -(prob_close[-1] - prob_wide[0]) / (wide_logP - close_logP)
                prob_intermediate = slope * (np.linspace(close_logP, wide_logP, neval) - close_logP) + prob_close[-1]
                prob_interp_int = interp1d(np.linspace(close_logP, wide_logP, neval), prob_intermediate)

                log_p_success = []
                n_success = 0
                while n_success < nsamp:
                    logP_samp = np.random.uniform(logP_lo_lim, logP_hi_lim, nsamp*5)
                    logP_prob = np.random.uniform(0, 1, nsamp*5)
        
                    logP_samp_lo = logP_samp[logP_samp<close_logP]
                    logP_prob_lo = logP_prob[logP_samp<close_logP]
                    log_p_success.extend(logP_samp_lo[np.where(logP_prob_lo < norm.pdf(logP_samp_lo, loc=mu, scale=sigma)*norm_close)])
        
                    logP_samp_int = logP_samp[(logP_samp>=close_logP) & (logP_samp<wide_logP)]
                    logP_prob_int = logP_prob[(logP_samp>=close_logP) & (logP_samp<wide_logP)]
                    log_p_success.extend(logP_samp_int[np.where(logP_prob_int < prob_interp_int(logP_samp_int))])
    
                    logP_samp_hi = logP_samp[(logP_samp>=wide_logP)]
                    logP_prob_hi = logP_prob[(logP_samp>=wide_logP)]

                    log_p_success.extend(logP_samp_hi[np.where(logP_prob_hi < norm.pdf(logP_samp_hi, loc=mu, scale=sigma)*norm_wide)])

                    n_success = len(log_p_success)
                log_p_success = np.array(log_p_success)[np.random.randint(0,n_success,nsamp)]    
                return log_p_success
            norm_wide, norm_close = utils.get_porb_norm(met)
            logP_dist = get_logP_dist(size, norm_wide, norm_close)
            logP_dist = logP_dist[np.random.randint(0, len(logP_dist), size)]
            porb = 10**logP_dist 
            aRL_over_a = a_min / utils.a_from_p(porb,mass1,mass2) 
            

        else:
            raise ValueError(
                "You have supplied a non-supported model; Please choose either log_flat, sana12, renzo19, raghavan10, or moe19"
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
            sana_max = np.array([0.9]*len(e_max))
            max_e = np.minimum(e_max, sana_max)
            ecc = utils.rndm(a=0.001, b=max_e, g=-0.45, size=size)
            
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

        from cosmic import _evolvebin


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
