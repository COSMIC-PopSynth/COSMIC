# -*- coding: utf-8 -*-
# Copyright (C) Katelyn Breivik (2017 - 2020)
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

from .sampler import register_sampler
from .independent import Sample
from .. import InitialBinaryTable

__author__ = 'Newlin Weatherford <newlinweatherford2017@u.northwestern.edu>'
__credits__ = ['Scott Coughlin <scott.coughlin@ligo.org>', 'Carl Rodriguez <carllouisrodriguez@gmail.com>']
__all__ = ['get_cmc_sampler', 'CMCSample']


def get_cmc_sampler(primary_model, ecc_model, porb_model, SF_start, SF_duration, binfrac_model, met, size, **kwargs):
    """Generates an initial binary sample according to user specified models

    Parameters
    ----------
    primary_model : `str`
        Model to sample primary mass; choices include: kroupa93, kroupa01, salpeter55

    ecc_model : `str`
        Model to sample eccentricity; choices include: thermal, uniform, sana12

    porb_model : `str`
        Model to sample orbital period; choices include: log_uniform, sana12

    SF_start : `float`
            Time in the past when star formation initiates in Myr

    SF_duration : `float`
            Duration of constant star formation beginning from SF_Start in Myr

    binfrac_model : `str or float`
        Model for binary fraction; choices include: vanHaaften or a fraction where 1.0 is 100% binaries

    met : `float`
        Sets the metallicity of the binary population where solar metallicity is 0.02

    size : `int`
        Size of the population to sample

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
    initconditions = CMCSample()

    # track the mass in singles and the mass in binaries
    mass_singles = 0.0
    mass_binaries = 0.0

    # track the total number of stars sampled
    n_singles = 0
    n_binaries = 0
    mass1, total_mass1 = initconditions.sample_primary(primary_model, size=size)
    mass1_binaries, mass_single, binfrac_binaries = initconditions.binary_select(mass1, binfrac_model=binfrac_model)
    mass2_binaries = initconditions.sample_secondary(mass1_binaries)

    # track the mass sampled
    mass_singles += np.sum(mass_single)
    mass_binaries += np.sum(mass1_binaries)
    mass_binaries += np.sum(mass2_binaries)

    # track the total number sampled
    n_singles += len(mass_single)
    n_binaries += len(mass1_binaries)

    # select out the primaries and secondaries that will produce the final kstars
    ind_select_primary, = np.where((mass1_binaries > primary_min) & (mass1_binaries < primary_max))
    ind_select_secondary, = np.where((mass2_binaries > secondary_min) & (mass2_binaries < secondary_max))
    ind_select = list(set(ind_select_primary).intersection(ind_select_secondary))
    mass1_binary.extend(mass1_binaries[ind_select])
    mass2_binary.extend(mass2_binaries[ind_select])
    binfrac.extend(binfrac_binaries[ind_select])

    mass1_binary = np.array(mass1_binary)
    mass2_binary = np.array(mass2_binary)
    binfrac = np.asarray(binfrac)
    ecc =  initconditions.sample_ecc(ecc_model, size = mass1_binary.size)
    porb =  initconditions.sample_porb(mass1_binary, mass2_binary, ecc, porb_model, size=mass1_binary.size)
    tphysf, metallicity = initconditions.sample_SFH(SF_start=SF_start, SF_duration=SF_duration, met=met, size = mass1_binary.size)
    metallicity[metallicity < 1e-4] = 1e-4
    metallicity[metallicity > 0.03] = 0.03
    kstar1 = initconditions.set_kstar(mass1_binary)
    kstar2 = initconditions.set_kstar(mass2_binary)

    return InitialBinaryTable.InitialCMCObjects(id_idx, initconditions.set_kstar(mass1), mass1, Reff, r, vr, vt, binind)


register_sampler('cmc', InitialBinaryTable, get_cmc_sampler,
                 usage="binfrac_model, primary_model, ecc_model, SFH_model, component_age, metallicity, size")


class CMCSample(Sample):
    def sample_vr(self):
        return
