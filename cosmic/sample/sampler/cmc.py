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

"""`cmc`
"""

import numpy as np

from .sampler import register_sampler
from .independent import Sample
from .. import InitialCMCTable, InitialBinaryTable
from ..cmc import elson
from ... import evolve

__author__ = 'Scott Coughlin <scott.coughlin@ligo.org>'
__credits__ = ['Newlin Weatherford <newlinweatherford2017@u.northwestern.edu>', 'Carl Rodriguez <carllouisrodriguez@gmail.com>']
__all__ = ['get_cmc_sampler', 'CMCSample']


def get_cmc_sampler(primary_model, ecc_model, porb_model, binfrac_model, met, size, **kwargs):
    """Generates an initial binary sample according to user specified models

    Parameters
    ----------
    primary_model : `str`
        Model to sample primary mass; choices include: kroupa93, kroupa01, salpeter55

    ecc_model : `str`
        Model to sample eccentricity; choices include: thermal, uniform, sana12

    porb_model : `str`
        Model to sample orbital period; choices include: log_uniform, sana12

    binfrac_model : `str or float`
        Model for binary fraction; choices include: vanHaaften or a fraction where 1.0 is 100% binaries

    met : `float`
        Sets the metallicity of the binary population where solar metallicity is 0.02

    size : `int`
        Size of the population to sample

    Returns
    -------
    InitialCMCTable : `pandas.DataFrame`
        DataFrame in the format of the InitialCMCTable

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
    mass1_binaries, mass_single, binfrac_binaries, binary_index = initconditions.binary_select(mass1, binfrac_model=binfrac_model)
    mass2_binaries = initconditions.sample_secondary(mass1_binaries)

    # track the mass sampled
    mass_singles += np.sum(mass_single)
    mass_binaries += np.sum(mass1_binaries)
    mass_binaries += np.sum(mass2_binaries)

    # track the total number sampled
    n_singles += len(mass_single)
    n_binaries += len(mass1_binaries)

    # select out the primaries and secondaries that will produce the final kstars
    mass1_binary = np.array(mass1_binaries)
    mass2_binary = np.array(mass2_binaries)
    binfrac = np.asarray(binfrac_binaries)
    ecc =  initconditions.sample_ecc(ecc_model, size = mass1_binary.size)
    porb =  initconditions.sample_porb(mass1_binary, mass2_binary, ecc, porb_model, size=mass1_binary.size)
    kstar1 = initconditions.set_kstar(mass1_binary)
    kstar2 = initconditions.set_kstar(mass2_binary)

    # obtain radius
    Reff = initconditions.set_reff(mass1, metallicity=met, **kwargs)
    Reff1 = Reff[binary_index]
    Reff2 = initconditions.set_reff(mass2_binary, metallicity=met, **kwargs)

    # set radial velocity, set transverse velocity, set location in cluster
    vr, vt, r = initconditions.set_vr_vt_r(N=mass1.size, **kwargs)

    # set singles id
    single_ids = np.arange(mass1.size)
    binary_secondary_object_id = np.arange(mass1.size, mass1.size + mass2_binaries.size)

    # set binind
    binind = np.zeros(mass1.size)
    binind[binary_index] = 1

    return InitialCMCTable.InitialCMCSingles(single_ids + 1, initconditions.set_kstar(mass1), mass1, Reff, r, vr, vt, binind), InitialCMCTable.InitialCMCBinaries(np.arange(mass1_binaries.size), single_ids[binary_index] + 1, kstar1, mass1_binary, Reff1, binary_secondary_object_id + 1, kstar2, mass2_binary, Reff2, porb, ecc)


register_sampler('cmc', InitialCMCTable, get_cmc_sampler,
                 usage="primary_model ecc_model porb_model binfrac_model met size")


class CMCSample(Sample):
    def set_vr_vt_r(self, **kwargs):
        cluster_profile = kwargs.pop('cluster_profile', 'elson')
        if cluster_profile == 'elson':
            elson_kwargs = {k:v for k,v in kwargs.items() if k in ['gamma', 'r_max', 'N']}
            vr, vt, r = elson.draw_vr_vt_r(**elson_kwargs)
        elif cluster_profile == 'plummer':
            plummer_kwargs = {k:v for k,v in kwargs.items() if k in ['r_max', 'N']}
            vr, vt, r = elson.draw_vr_vt_r(gamma=4, **plummer_kwargs)
        else:
            raise ValueError("Cluster profile passed not defined")

        return vr, vt, r

    def set_reff(self, mass, metallicity, **kwargs):
        # NUMBER 1: PASS A DICTIONARY OF FLAGS
        BSEDict = kwargs.pop('BSEDict', {})

        # NUMBER 2: PASS PATH TO A INI FILE WITH THE FLAGS DEFINED
        params = kwargs.pop('params', None)

        nproc = min(kwargs.pop('nproc', 1), mass.size)

        bpp, bcm, initial_conditions, kick_info = evolve.Evolve.evolve(InitialBinaryTable.InitialBinaries(mass, np.ones_like(mass)*0, np.ones_like(mass)*-1, np.ones_like(mass)*-1, np.ones_like(mass)*0.1, self.set_kstar(mass), np.ones_like(mass)*0, np.ones_like(mass)*metallicity), BSEDict=BSEDict, params=params, nproc=nproc)

        return bcm.groupby('bin_num').first().rad_1
