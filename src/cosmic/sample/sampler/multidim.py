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

"""`multidim`
"""
from schwimmbad import MultiPool, MPIPool

from .sampler import register_sampler
from .. import InitialBinaryTable
from ... import utils

import numpy as np
import pandas as pd

__author__ = "Katelyn Breivik <katie.breivik@gmail.com>"
__credits__ = "Scott Coughlin <scott.coughlin@ligo.org>"
__all__ = ["get_multidim_sampler", "MultiDim"]


def get_multidim_sampler(
    final_kstar1,
    final_kstar2,
    rand_seed,
    nproc,
    SF_start,
    SF_duration,
    met,
    size,
    **kwargs
):
    """adapted version of Maxwell Moe's IDL code that generates a population of single and binary stars

    Below is the adapted version of Maxwell Moe's IDL code
    that generates a population of single and binary stars
    based on the paper Mind your P's and Q's
    By Maxwell Moe and Rosanne Di Stefano

    The python code has been adopted by Mads Sørensen

    Version history:
    V. 0.1; 2017/02/03
    By Mads Sørensen
    - This is a pure adaption from IDL to Python.
    - The function idl_tabulate is similar to
    the IDL function int_tabulated except, this function seems to be slightly
    more exact in its solution.
    Therefore, relative to the IDL code, there are small numerical differences.

    Comments below beginning with ; is the original nodes by Maxwell Moe.
    Please read these careful for understanding the script.
    ; NOTE - This version produces only the statistical distributions of
    ;        single stars, binaries, and inner binaries in hierarchical triples.
    ;        Outer tertiaries in hierarchical triples are NOT generated.
    ;        Moreover, given a set of companions, all with period P to
    ;        primary mass M1, this version currently uses an approximation to
    ;        determine the fraction of those companions that are inner binaries
    ;        vs. outer triples. Nevertheless, this approximation reproduces
    ;        the overall multiplicity statistics.
    ; Step 1 - Tabulate probably density functions of periods,
    ;          mass ratios, and eccentricities based on
    ;          analytic fits to corrected binary star populations.
    ; Step 2 - Implement Monte Carlo method to generate stellar
    ;          population from those density functions.

    Parameters
    ----------
    final_kstar1 : `list` or `int`
        Int or list of final kstar1

    final_kstar2 : `list` or `int`
        Int or list of final kstar2

    rand_seed : `int`
        Int to seed random number generator

    nproc : `int`
        Number of processors to use to generate population

    SF_start : `float`
            Time in the past when star formation initiates in Myr

    SF_duration : `float`
            Duration of constant star formation beginning from SF_Start in Myr

    met : `float`
        Sets the metallicity of the binary population where solar metallicity is 0.02

    size : `int`
        Size of the population to sample

    **porb_lo : `float`
        Lower limit in days for the orbital period distribution

    **porb_hi: `float`
        Upper limit in days for the orbital period distribution

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
    porb_lo = kwargs.pop("porb_lo", 0.15)
    porb_hi = kwargs.pop("porb_hi", 8.0)
    pool = kwargs.pop("pool", None)
    mp_seeds = kwargs.pop("mp_seeds", None)

    primary_min, primary_max, secondary_min, secondary_max = utils.mass_min_max_select(
        final_kstar1, final_kstar2
    )

    initconditions = MultiDim()

    (
        mass1_binary,
        mass2_binary,
        porb,
        ecc,
        single_mass_list,
        mass_singles,
        mass_binaries,
        n_singles,
        n_binaries,
        binfrac,
    ) = initconditions.initial_sample(
        primary_min,
        secondary_min,
        primary_max,
        secondary_max,
        porb_lo,
        porb_hi,
        rand_seed,
        size=size,
        nproc=nproc,
        pool=pool,
        mp_seeds=mp_seeds,
    )

    tphysf, metallicity = initconditions.sample_SFH(
        SF_start=SF_start, SF_duration=SF_duration, met=met, size=mass1_binary.size
    )

    kstar1 = initconditions.set_kstar(mass1_binary)
    kstar2 = initconditions.set_kstar(mass2_binary)
    metallicity[metallicity < 1e-4] = 1e-4
    metallicity[metallicity > 0.03] = 0.03

    if kwargs.pop("keep_singles", True):
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
            SF_start=SF_start, SF_duration=SF_duration, met=met, size=single_mass_list.size
        )
        metallicity[metallicity < 1e-4] = 1e-4
        metallicity[metallicity > 0.03] = 0.03
        kstar1 = initconditions.set_kstar(single_mass_list)
        singles_table = InitialBinaryTable.InitialBinaries(
            single_mass_list,
            np.ones_like(single_mass_list)*0,
            np.ones_like(single_mass_list)*-1,
            np.ones_like(single_mass_list)*-1,
            tphysf,
            kstar1,
            np.ones_like(single_mass_list)*15, # # kstar2 is not used for singles
            metallicity,
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
        mass_singles,
        mass_binaries,
        n_singles,
        n_binaries,
    )

register_sampler(
    "multidim",
    InitialBinaryTable,
    get_multidim_sampler,
    usage="final_kstar1, final_kstar2, rand_seed, nproc, SFH_model, component_age, metallicity, size, binfrac",
)


class MultiDim:

    # -----------------------------------
    # Below is the adapted version of Maxwell Moe's IDL code
    # that generates a population of single and binary stars
    # based on the paper Mind your P's and Q's
    # By Maxwell Moe and Rosanne Di Stefano
    #
    # The python code has been adopted by Mads Sørensen
    # -----------------------------------
    # Version history:
    # V. 0.1; 2017/02/03
    # By Mads Sørensen
    # - This is a pure adaption from IDL to Python.
    # - The function idl_tabulate is similar to
    # the IDL function int_tabulated except, this function seems to be slightly
    # more exact in its solution.
    # Therefore, relative to the IDL code, there are small numerical differences.
    # -----------------------------------

    #
    # Comments below beginning with ; is the original nodes by Maxwell Moe.
    # Please read these careful for understanding the script.
    # ; NOTE - This version produces only the statistical distributions of
    # ;        single stars, binaries, and inner binaries in hierarchical triples.
    # ;        Outer tertiaries in hierarchical triples are NOT generated.
    # ;        Moreover, given a set of companions, all with period P to
    # ;        primary mass M1, this version currently uses an approximation to
    # ;        determine the fraction of those companions that are inner binaries
    # ;        vs. outer triples. Nevertheless, this approximation reproduces
    # ;        the overall multiplicity statistics.
    # ; Step 1 - Tabulate probably density functions of periods,
    # ;          mass ratios, and eccentricities based on
    # ;          analytic fits to corrected binary star populations.
    # ; Step 2 - Implement Monte Carlo method to generate stellar
    # ;          population from those density functions.
    # ;
    def initial_sample(
        self,
        M1min=0.08,
        M2min=0.08,
        M1max=150.0,
        M2max=150.0,
        porb_lo=0.15,
        porb_hi=8.0,
        rand_seed=0,
        size=None,
        nproc=1,
        pool=None,
        mp_seeds=None,
    ):
        """Sample initial binary distribution according to Moe & Di Stefano (2017)
        <http://adsabs.harvard.edu/abs/2017ApJS..230...15M>`_

        Parameters
        ----------
        M1min : `float`
            minimum primary mass to sample [Msun]
            DEFAULT: 0.08
        M2min : `float`
            minimum secondary mass to sample [Msun]
            DEFAULT: 0.08
        M1max : `float`
            maximum primary mass to sample [Msun]
            DEFAULT: 150.0
        M2max : `float`
            maximum primary mass to sample [Msun]
            DEFAULT: 150.0
        porb_lo : `float`
            minimum orbital period to sample [log10(days)]
        porb_hi : `float`
            maximum orbital period to sample [log10(days)]
        rand_seed : int
            random seed generator
            DEFAULT: 0
        size : int, optional
            number of evolution times to sample
            NOTE: this is set in cosmic-pop call as Nstep

        Returns
        -------
        primary_mass_list : array
            array of primary masses with size=size
        secondary_mass_list : array
            array of secondary masses with size=size
        porb_list : array
            array of orbital periods in days with size=size
        ecc_list : array
            array of eccentricities with size=size
        single_mass_list : array
            array of mass of single stars
        mass_singles : `float`
            Total mass in single stars needed to generate population
        mass_binaries : `float`
            Total mass in binaries needed to generate population
        n_singles : `int`
            Number of single stars needed to generate a population
        n_binaries : `int`
            Number of binaries needed to generate a population
        binfrac_list : array
            array of binary probabilities based on primary mass and period with size=size
        """
        if pool is None:
            with MultiPool(processes=nproc) as pool:
                if mp_seeds is not None:
                    if len(list(mp_seeds)) != nproc:
                        raise ValueError("Must supply a list of random seeds with length equal to number of processors")
                else:
                    mp_seeds = [nproc * (task._identity[0]-1) for task in pool._pool]

                inputs = [(M1min, M2min, M1max, M2max, porb_hi, porb_lo, size/nproc, rand_seed + mp_seed)
                          for mp_seed in mp_seeds]
                worker = Worker()
                results = list(pool.map(worker, inputs))
        else:
            if mp_seeds is not None:
                if len(list(mp_seeds)) != nproc:
                    raise ValueError("Must supply a list of random seeds with length equal to number of processors")
            else:
                if isinstance(pool, MPIPool):
                    mp_seeds = [nproc * (task - 1) for task in pool.workers]
                elif isinstance(pool, MultiPool):
                    mp_seeds = [nproc * (task._identity[0] - 1) for task in pool._pool]
                else:
                    mp_seeds = [0 for i in range(nproc)]

            inputs = [(M1min, M2min, M1max, M2max, porb_hi, porb_lo, size/nproc, rand_seed + mp_seed) for mp_seed in mp_seeds]
            worker = Worker()
            results = list(pool.map(worker, inputs))

        dat_lists = [[], [], [], [], [], [], [], [], [], []]

        for output_list in results:
            ii = 0
            for dat_list in output_list:
                dat_lists[ii].append(dat_list)
                ii += 1

        primary_mass_list = np.hstack(dat_lists[0])
        secondary_mass_list = np.hstack(dat_lists[1])
        porb_list = np.hstack(dat_lists[2])
        ecc_list = np.hstack(dat_lists[3])
        single_mass_list = np.hstack(dat_lists[4])
        mass_singles = np.sum(dat_lists[5])
        mass_binaries = np.sum(dat_lists[6])
        n_singles = np.sum(dat_lists[7])
        n_binaries = np.sum(dat_lists[8])
        binfrac_list = np.hstack(dat_lists[9])

        return (
                primary_mass_list,
                secondary_mass_list,
                porb_list,
                ecc_list,
                single_mass_list,
                mass_singles,
                mass_binaries,
                n_singles,
                n_binaries,
                binfrac_list
               )

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
            metallicity = np.ones(size)*met
            return tphys, metallicity
        else:
            raise ValueError('SF_start and SF_duration must be positive and SF_start must be greater than 0.0')

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


class Worker(object):
    def __call__(self, task):
        M1min, M2min, M1max, M2max, porb_hi, porb_lo, size, seed = task
        return self._sample_initial_pop(M1min, M2min, M1max, M2max, porb_hi, porb_lo, size, seed)

    def _sample_initial_pop(self, M1min, M2min, M1max, M2max, porb_hi, porb_lo, size, seed):
        # Tabulate probably density functions of periods,
        # mass ratios, and eccentricities based on
        # analytic fits to corrected binary star populations.

        numM1 = 101
        # use binwidths to maintain structure of original array
        # default size is: numlogP=158
        bwlogP = 0.05
        numq = 91
        nume = 100

        # ; Vector of primary masses M1 (Msun), logarithmic orbital period P (days),
        # ; mass ratios q = Mcomp/M1, and eccentricities e
        #
        # ; 0.8 < M1 < 40 (where we have statistics corrected for selection effects)
        M1_lo = 0.8
        M1_hi = 40

        M1v = np.logspace(np.log10(M1_lo), np.log10(M1_hi), numM1)
        # ; 0.15 < log P < 8.0
        # ; or use user specified values
        log10_porb_lo = porb_lo
        log10_porb_hi = porb_hi
        logPv = np.arange(log10_porb_lo, log10_porb_hi + bwlogP, bwlogP)
        numlogP = len(logPv)

        # ; 0.10 < q < 1.00
        q_lo = 0.1
        q_hi = 1.0
        qv = np.linspace(q_lo, q_hi, numq)

        # ; 0.0001 < e < 0.9901
        # ; set minimum to non-zero value to avoid numerical errors
        e_lo = 0.0
        e_hi = 0.99
        ev = np.linspace(e_lo, e_hi, nume) + 0.0001
        # ; Note that companions outside this parameter space (e.g., q < 0.1,
        # ; log P (days) > 8.0 are not constrained in M+D16 and therefore
        # ; not considered.

        # ; Distribution functions - define here, but evaluate within for loops.

        # ; Frequency of companions with q > 0.1 per decade of orbital period.
        # ; Bottom panel in Fig. 37 of M+D17
        flogP_sq = np.zeros([numlogP, numM1])

        # ; Given M1 and P, the cumulative distribution of mass ratios q
        cumqdist = np.zeros([numq, numlogP, numM1])

        # ; Given M1 and P, the cumulative distribution of eccentricities e
        cumedist = np.zeros([nume, numlogP, numM1])

        # ; Given M1 and P, the probability that the companion
        # ; is a member of the inner binary (currently an approximation).
        # ; 100% for log P < 1.5, decreases with increasing P
        probbin = np.zeros([numlogP, numM1])

        # ; Given M1, the cumulative period distribution of the inner binary
        # ; Normalized so that max(cumPbindist) = total binary frac. (NOT unity)
        cumPbindist = np.zeros([numlogP, numM1])
        # ; Slope alpha of period distribution across intermediate periods
        # ; 2.7 - DlogP < log P < 2.7 + DlogP, see Section 9.3 and Eqn. 23.
        # ; Slightly updated from version 1.
        alpha = 0.018
        DlogP = 0.7

        # ; Heaviside function for twins with 0.95 < q < 1.00
        H = np.zeros(numq)
        ind = np.where(qv >= 0.95)
        H[ind] = 1.0
        H = H / utils.idl_tabulate(qv, H)  # ;normalize so that integral is unity

        # ; Relevant indices with respect to mass ratio
        indlq = np.where(qv >= 0.3)
        indsq = np.where(qv < 0.3)
        indq0p3 = np.min(indlq)

        # FILL IN THE MULTIDIMENSIONAL DISTRIBUTION FUNCTIONS
        # ; Loop through primary mass
        for i in range(0, numM1):
            myM1 = M1v[i]
            # ; Twin fraction parameters that are dependent on M1 only; section 9.1
            FtwinlogPle1 = 0.3 - 0.15 * np.log10(myM1)  # ; Eqn. 6
            logPtwin = 8.0 - myM1  # ; Eqn. 7a
            if myM1 >= 6.5:
                logPtwin = 1.5  # ; Eqn. 7b
            # ; Frequency of companions with q > 0.3 at different orbital periods
            # ; and dependent on M1 only; section 9.3 (slightly modified since v1)
            flogPle1 = (
                0.020 + 0.04 * np.log10(myM1) + 0.07 * (np.log10(myM1)) ** 2.0
            )  # ; Eqn. 20
            flogPeq2p7 = (
                0.039 + 0.07 * np.log10(myM1) + 0.01 * (np.log10(myM1)) ** 2.0
            )  # ; Eqn. 21
            flogPeq5p5 = (
                0.078 - 0.05 * np.log10(myM1) + 0.04 * (np.log10(myM1)) ** 2.0
            )  # ; Eqn. 22
            # ; Loop through orbital period P
            for j in range(0, numlogP):
                mylogP = logPv[j]
                # ; Given M1 and P, set excess twin fraction; section 9.1 and Eqn. 5
                if mylogP <= 1.0:
                    Ftwin = FtwinlogPle1
                else:
                    Ftwin = FtwinlogPle1 * (1.0 - (mylogP - 1.0) / (logPtwin - 1.0))
                if mylogP >= logPtwin:
                    Ftwin = 0.0

                # ; Power-law slope gamma_largeq for M1 < 1.2 Msun and various P; Eqn. 9
                if mylogP <= 5.0:
                    gl_1p2 = -0.5
                if mylogP > 5.0:
                    gl_1p2 = -0.5 - 0.3 * (mylogP - 5.0)

                # ; Power-law slope gamma_largeq for M1 = 3.5 Msun and various P; Eqn. 10
                if mylogP <= 1.0:
                    gl_3p5 = -0.5
                if (mylogP > 1.0) and (mylogP <= 4.5):
                    gl_3p5 = -0.5 - 0.2 * (mylogP - 1.0)
                if (mylogP > 4.5) and (mylogP <= 6.5):
                    gl_3p5 = -1.2 - 0.4 * (mylogP - 4.5)
                if mylogP > 6.5:
                    gl_3p5 = -2.0

                # ; Power-law slope gamma_largeq for M1 > 6 Msun and various P; Eqn. 11
                if mylogP <= 1.0:
                    gl_6 = -0.5
                if (mylogP > 1.0) and (mylogP <= 2.0):
                    gl_6 = -0.5 - 0.9 * (mylogP - 1.0)
                if (mylogP > 2.0) and (mylogP <= 4.0):
                    gl_6 = -1.4 - 0.3 * (mylogP - 2.0)

                if mylogP > 4.0:
                    gl_6 = -2.0

                # ; Given P, interpolate gamma_largeq w/ respect to M1 at myM1
                if myM1 <= 1.2:
                    gl = gl_1p2
                if (myM1 > 1.2) and (myM1 <= 3.5):
                    gl = np.interp(
                        np.log10(myM1), np.log10([1.2, 3.5]), [gl_1p2, gl_3p5]
                    )
                if (myM1 > 3.5) and (myM1 <= 6.0):
                    gl = np.interp(np.log10(myM1), np.log10([3.5, 6.0]), [gl_3p5, gl_6])
                if myM1 > 6.0:
                    gl = gl_6

                # ; Power-law slope gamma_smallq for M1 < 1.2 Msun and all P; Eqn. 13
                gs_1p2 = 0.3

                # ; Power-law slope gamma_smallq for M1 = 3.5 Msun and various P; Eqn. 14
                if mylogP <= 2.5:
                    gs_3p5 = 0.2
                if (mylogP > 2.5) and (mylogP <= 5.5):
                    gs_3p5 = 0.2 - 0.3 * (mylogP - 2.5)
                if mylogP > 5.5:
                    gs_3p5 = -0.7 - 0.2 * (mylogP - 5.5)

                # ; Power-law slope gamma_smallq for M1 > 6 Msun and various P; Eqn. 15
                if mylogP <= 1.0:
                    gs_6 = 0.1
                if (mylogP > 1.0) and (mylogP <= 3.0):
                    gs_6 = 0.1 - 0.15 * (mylogP - 1.0)
                if (mylogP > 3.0) and (mylogP <= 5.6):
                    gs_6 = -0.2 - 0.50 * (mylogP - 3.0)
                if mylogP > 5.6:
                    gs_6 = -1.5

                # ; Given P, interpolate gamma_smallq w/ respect to M1 at myM1
                if myM1 <= 1.2:
                    gs = gs_1p2
                if (myM1 > 1.2) and (myM1 <= 3.5):
                    gs = np.interp(
                        np.log10(myM1), np.log10([1.2, 3.5]), [gs_1p2, gs_3p5]
                    )
                if (myM1 > 3.5) and (myM1 <= 6.0):
                    gs = np.interp(np.log10(myM1), np.log10([3.5, 6.0]), [gs_3p5, gs_6])
                if myM1 > 6.0:
                    gs = gs_6

                # ; Given Ftwin, gamma_smallq, and gamma_largeq at the specified M1 & P,
                # ; tabulate the cumulative mass ratio distribution across 0.1 < q < 1.0
                fq = qv ** gl  # ; slope across 0.3 < q < 1.0
                fq = fq / utils.idl_tabulate(
                    qv[indlq], fq[indlq]
                )  # ; normalize to 0.3 < q < 1.0
                fq = fq * (1.0 - Ftwin) + H * Ftwin  # ; add twins
                fq[indsq] = (
                    fq[indq0p3] * (qv[indsq] / 0.3) ** gs
                )  # ; slope across 0.1 < q < 0.3
                cumfq = np.cumsum(fq) - fq[0]  # ; cumulative distribution
                cumfq = cumfq / np.max(cumfq)  # ; normalize cumfq(q=1.0) = 1
                cumqdist[:, j, i] = cumfq  # ; save to grid

                # ; Given M1 and P, q_factor is the ratio of all binaries 0.1 < q < 1.0
                # ; to those with 0.3 < q < 1.0
                q_factor = utils.idl_tabulate(qv, fq)

                # ; Given M1 & P, calculate power-law slope eta of eccentricity dist.
                if mylogP >= 0.7:
                    # ; For log P > 0.7 use fits in Section 9.2.
                    # ; Power-law slope eta for M1 < 3 Msun and log P > 0.7
                    eta_3 = 0.6 - 0.7 / (mylogP - 0.5)  # ; Eqn. 17
                    # ; Power-law slope eta for M1 > 7 Msun and log P > 0.7
                    eta_7 = 0.9 - 0.2 / (mylogP - 0.5)  # ; Eqn. 18
                else:
                    # ; For log P < 0.7, set eta to fitted values at log P = 0.7
                    eta_3 = -2.9
                    eta_7 = -0.1

                # ; Given P, interpolate eta with respect to M1 at myM1
                if myM1 <= 3.0:
                    eta = eta_3
                if (myM1 > 3.0) and (myM1 <= 7.0):
                    eta = np.interp(
                        np.log10(myM1), np.log10([3.0, 7.0]), [eta_3, eta_7]
                    )
                if myM1 > 7.0:
                    eta = eta_7

                # ; Given eta at the specified M1 and P, tabulate eccentricity distribution
                if 10 ** mylogP <= 2.0:
                    # ; For P < 2 days, assume all systems are close to circular
                    # ; For adopted ev (spacing and minimum value), eta = -3.2 satisfies this
                    fe = ev ** (-3.2)
                else:
                    fe = ev ** eta
                    e_max = 1.0 - (10 ** mylogP / 2.0) ** (
                        -2.0 / 3.0
                    )  # ; maximum eccentricity for given P
                    ind = np.where(ev >= e_max)
                    fe[ind] = 0.0  # ; set dist. = 0 for e > e_max
                    # ; Assume e dist. has power-law slope eta for 0.0 < e / e_max < 0.8 and
                    # ; then linear turnover between 0.8 < e / e_max < 1.0 so that dist.
                    # ; is continuous at e / e_max = 0.8 and zero at e = e_max
                    ind = np.where((ev >= 0.8 * e_max) & (ev <= 1.0 * e_max))
                    ind_cont = np.min(ind) - 1
                    fe[ind] = np.interp(
                        ev[ind], [0.8 * e_max, 1.0 * e_max], [fe[ind_cont], 0.0]
                    )

                cumfe = np.cumsum(fe) - fe[0]  # ; cumulative distribution
                cumfe = cumfe / np.max(cumfe)  # ; normalize cumfe(e=e_max) = 1
                cumedist[:, j, i] = cumfe  # ; save to grid

                # ; Given constants alpha and DlogP and
                # ; M1 dependent values flogPle1, flogPeq2p7, and flogPeq5p5,
                # ; calculate frequency flogP of companions with q > 0.3 per decade
                # ; of orbital period at given P (Section 9.3 and Eqn. 23)
                if mylogP <= 1.0:
                    flogP = flogPle1
                if (mylogP > 1.0) and (mylogP <= 2.7 - DlogP):
                    flogP = flogPle1 + (mylogP - 1.0) / (1.7 - DlogP) * (
                        flogPeq2p7 - flogPle1 - alpha * DlogP
                    )
                if (mylogP > 2.7 - DlogP) and (mylogP <= 2.7 + DlogP):
                    flogP = flogPeq2p7 + alpha * (mylogP - 2.7)
                if (mylogP > 2.7 + DlogP) and (mylogP <= 5.5):
                    flogP = (
                        flogPeq2p7
                        + alpha * DlogP
                        + (mylogP - 2.7 - DlogP)
                        / (2.8 - DlogP)
                        * (flogPeq5p5 - flogPeq2p7 - alpha * DlogP)
                    )
                if mylogP > 5.5:
                    flogP = flogPeq5p5 * np.exp(-0.3 * (mylogP - 5.5))

                # ; Convert frequency of companions with q > 0.3 to frequency of
                # ; companions with q > 0.1 according to q_factor; save to grid
                flogP_sq[j, i] = flogP * q_factor

                # ; Calculate prob. that a companion to M1 with period P is the
                # ; inner binary.  Currently this is an approximation.
                # ; 100% for log P < 1.5
                # ; For log P > 1.5 adopt functional form that reproduces M1 dependent
                # ; multiplicity statistics in Section 9.4, including a
                # ; 41% binary star faction (59% single star fraction) for M1 = 1 Msun and
                # ; 96% binary star fraction (4% single star fraction) for M1 = 28 Msun
                if mylogP <= 1.5:
                    probbin[j, i] = 1.0
                else:
                    probbin[j, i] = (
                        1.0 - 0.11 * (mylogP - 1.5) ** 1.43 * (myM1 / 10.0) ** 0.56
                    )
                if probbin[j, i] <= 0.0:
                    probbin[j, i] = 0.0

            # ; Given M1, calculate cumulative binary period distribution
            mycumPbindist = (
                np.cumsum(flogP_sq[:, i] * probbin[:, i])
                - flogP_sq[0, i] * probbin[0, i]
            )
            # ; Normalize so that max(cumPbindist) = total binary star fraction (NOT 1)
            mycumPbindist = (
                mycumPbindist
                / np.max(mycumPbindist)
                * utils.idl_tabulate(logPv, flogP_sq[:, i] * probbin[:, i])
            )
            cumPbindist[:, i] = mycumPbindist  # ;save to grid

        # ; Step 2
        # ; Implement Monte Carlo method / random number generator to select
        # ; single stars and binaries from the grids of distributions

        # ; Create vector for PRIMARY mass function, which is the mass distribution
        # ; of single stars and primaries in binaries.
        # ; This is NOT the IMF, which is the mass distribution of single stars,
        # ; primaries in binaries, and secondaries in binaries.

        np.random.seed(seed)

        mass_singles = 0.0
        mass_binaries = 0.0
        n_singles = 0
        n_binaries = 0
        primary_mass_list = []
        secondary_mass_list = []
        single_mass_list = []
        porb_list = []
        ecc_list = []
        binfrac_list = []

        # Full primary mass vector across 0.08 < M1 < 150
        M1 = np.linspace(0, 150, 150000) + 0.08
        # Slope = -2.3 for M1 > 1 Msun
        fM1 = M1**(-2.3)
        # Slope = -1.6 for M1 = 0.5 - 1.0 Msun
        ind = np.where(M1 <= 1.0)
        fM1[ind] = M1[ind]**(-1.6)
        # Slope = -0.8 for M1 = 0.15 - 0.5 Msun
        ind = np.where(M1 <= 0.5)
        fM1[ind] = M1[ind]**(-0.8) / 0.5**(1.6 - 0.8)
        # Cumulative primary mass distribution function
        cumfM1 = np.cumsum(fM1) - fM1[0]
        cumfM1 = cumfM1 / np.max(cumfM1)
        # Value of primary mass CDF where M1 = M1min
        # Minimum primary mass to generate (must be >0.080 Msun)
        cumf_M1min = np.interp(0.08, M1, cumfM1)
        while len(primary_mass_list) < size:

            # Select primary M1 > M1min from primary mass function
            myM1 = np.interp(cumf_M1min + (1.0 - cumf_M1min) * np.random.rand(), cumfM1, M1)

            # Find index of M1v that is closest to myM1.
            #     ; For M1 = 40 - 150 Msun, adopt binary statistics of M1 = 40 Msun.
            #     ; For M1 = 0.08 - 0.8 Msun, adopt P and e dist of M1 = 0.8Msun,
            #     ; scale and interpolate the companion frequencies so that the
            #     ; binary star fraction of M1 = 0.08 Msun primaries is zero,
            #     ; and truncate the q distribution so that q > q_min = 0.08/M1
            indM1 = np.where(abs(myM1 - M1v) == min(abs(myM1 - M1v)))
            indM1 = indM1[0]

            # ; Given M1, determine cumulative binary period distribution
            mycumPbindist_flat = (cumPbindist[:, indM1]).flatten()
            # If M1 < 0.8 Msun, rescale to appropriate binary star fraction
            if(myM1 <= 0.8):
                mycumPbindist_flat = mycumPbindist_flat * np.interp(np.log10(myM1), np.log10([0.08, 0.8]), [0.0, 1.0])

            # ; Given M1, determine the binary star fraction
            mybinfrac = np.max(mycumPbindist_flat)

            # ; Generate random number myrand between 0 and 1
            myrand = np.random.rand()
            # If random number < binary star fraction, generate a binary
            if(myrand < mybinfrac):
                # Given myrand, select P and corresponding index in logPv
                mylogP = np.interp(myrand, mycumPbindist_flat, logPv)
                indlogP = np.where(abs(mylogP - logPv) == min(abs(mylogP - logPv)))
                indlogP = indlogP[0]

                # Given M1 & P, select e from eccentricity distribution
                mye = np.interp(np.random.rand(), cumedist[:, indlogP, indM1].flatten(), ev)

                # Given M1 & P, determine mass ratio distribution.
                # If M1 < 0.8 Msun, truncate q distribution and consider
                # only mass ratios q > q_min = 0.08 / M1
                mycumqdist = cumqdist[:, indlogP, indM1].flatten()
                if(myM1 < 0.8):
                    q_min = 0.08 / myM1
                    # Calculate cumulative probability at q = q_min
                    cum_qmin = np.interp(q_min, qv, mycumqdist)
                    # Rescale and renormalize cumulative distribution for q > q_min
                    mycumqdist = mycumqdist - cum_qmin
                    mycumqdist = mycumqdist / max(mycumqdist)
                    # Set probability = 0 where q < q_min
                    indq = np.where(qv <= q_min)
                    mycumqdist[indq] = 0.0

                # Given M1 & P, select q from cumulative mass ratio distribution
                myq = np.interp(np.random.rand(), mycumqdist, qv)

                if ((myM1 > M1min) and (myq * myM1 > M2min) and (myM1 < M1max) and
                   (myq * myM1 < M2max) and (mylogP < porb_hi) and (mylogP > porb_lo)):
                    primary_mass_list.append(myM1)
                    secondary_mass_list.append(myq * myM1)
                    porb_list.append(10**mylogP)
                    ecc_list.append(mye)
                    binfrac_list.append(mybinfrac)
                mass_binaries += myM1
                mass_binaries += myq * myM1
                n_binaries += 1
            else:
                single_mass_list.append(myM1)
                mass_singles += myM1
                n_singles += 1

        return (
                primary_mass_list,
                secondary_mass_list,
                porb_list,
                ecc_list,
                single_mass_list,
                mass_singles,
                mass_binaries,
                n_singles,
                n_binaries,
                binfrac_list
               )
