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
from ... import _evolvebin
from ... import utils

__author__ = "Scott Coughlin <scott.coughlin@ligo.org>"
__credits__ = [
    "Newlin Weatherford <newlinweatherford2017@u.northwestern.edu>",
    "Carl Rodriguez <carllouisrodriguez@gmail.com>",
]
__all__ = ["get_cmc_sampler", "CMCSample"]


def get_cmc_sampler(
    primary_model, ecc_model, porb_model, qmin, binfrac_model, met, size, **kwargs
):
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

    qmin : `float`
        Sets the minimum mass ratio if q>0 and uses the pre-MS lifetime
        of the companion for primary mass > 5 msun to set q and sets q=0.1 
        for all primary mass < 5 msun if q < 0

        cluster_profile : `str`
        Model to use for the cluster profile (i.e. sampling of the placement of objects in the cluster and their velocity within the cluster)
        options include king, plummer and elson.

    params : `str`
        Path to the inifile with the BSE parameters. We need to generate radii for the single stars of the cluster by
        running BSE for a tiny time step.

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
    (
        mass1_binaries,
        mass_single,
        binfrac_binaries,
        binary_index,
    ) = initconditions.binary_select(mass1, binfrac_model=binfrac_model)
    mass2_binaries = initconditions.sample_secondary(mass1_binaries, qmin)

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
    ecc = initconditions.sample_ecc(ecc_model, size=mass1_binary.size)
    porb = initconditions.sample_porb(
        mass1_binary, mass2_binary, ecc, porb_model, size=mass1_binary.size
    )
    sep = utils.a_from_p(porb, mass1_binary, mass2_binary)
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

    # set binind and correct masses of binaries
    binind = np.zeros(mass1.size)
    binind[binary_index] = np.arange(len(binary_index)) + 1
    mass1[binary_index] += mass2_binary

    singles_table = InitialCMCTable.InitialCMCSingles(
        single_ids + 1, initconditions.set_kstar(mass1), mass1, Reff, r, vr, vt, binind
    )
    singles_table.metallicity = met
    singles_table.mass_of_cluster = np.sum(singles_table["m"])

    binaries_table = InitialCMCTable.InitialCMCBinaries(
        np.arange(mass1_binaries.size) + 1,
        single_ids[binary_index] + 1,
        kstar1,
        mass1_binary,
        Reff1,
        binary_secondary_object_id + 1,
        kstar2,
        mass2_binary,
        Reff2,
        sep,
        ecc,
    )
    binaries_table.metallicity = met

    return singles_table, binaries_table


register_sampler(
    "cmc",
    InitialCMCTable,
    get_cmc_sampler,
    usage="primary_model ecc_model porb_model binfrac_model met size",
)


class CMCSample(Sample):
    def set_vr_vt_r(self, **kwargs):
        cluster_profile = kwargs.pop("cluster_profile", "elson")
        if cluster_profile == "elson":
            elson_kwargs = {
                k: v for k, v in kwargs.items() if k in ["gamma", "r_max", "N"]
            }
            vr, vt, r = elson.draw_vr_vt_r(**elson_kwargs)
        elif cluster_profile == "plummer":
            plummer_kwargs = {k: v for k, v in kwargs.items() if k in ["r_max", "N"]}
            vr, vt, r = elson.draw_vr_vt_r(gamma=4, **plummer_kwargs)
        else:
            raise ValueError("Cluster profile passed not defined")

        return vr, vt, r

    def set_reff(self, mass, metallicity, **kwargs):
        # NUMBER 1: PASS A DICTIONARY OF FLAGS
        BSEDict = kwargs.pop("BSEDict", {})

        # NUMBER 2: PASS PATH TO A INI FILE WITH THE FLAGS DEFINED
        params = kwargs.pop("params", None)

        if params is not None:
            BSEDict, _, _, _, _ = utils.parse_inifile(params)

        # set BSE consts
        _evolvebin.windvars.neta = BSEDict["neta"]
        _evolvebin.windvars.bwind = BSEDict["bwind"]
        _evolvebin.windvars.hewind = BSEDict["hewind"]
        _evolvebin.cevars.alpha1 = BSEDict["alpha1"]
        _evolvebin.cevars.lambdaf = BSEDict["lambdaf"]
        _evolvebin.ceflags.ceflag = BSEDict["ceflag"]
        _evolvebin.flags.tflag = BSEDict["tflag"]
        _evolvebin.flags.ifflag = BSEDict["ifflag"]
        _evolvebin.flags.wdflag = BSEDict["wdflag"]
        _evolvebin.snvars.pisn = BSEDict["pisn"]
        _evolvebin.flags.bhflag = BSEDict["bhflag"]
        _evolvebin.flags.remnantflag = BSEDict["remnantflag"]
        _evolvebin.ceflags.cekickflag = BSEDict["cekickflag"]
        _evolvebin.ceflags.cemergeflag = BSEDict["cemergeflag"]
        _evolvebin.ceflags.cehestarflag = BSEDict["cehestarflag"]
        _evolvebin.flags.grflag = BSEDict["grflag"]
        _evolvebin.flags.bhms_coll_flag = BSEDict["bhms_coll_flag"]
        _evolvebin.snvars.mxns = BSEDict["mxns"]
        _evolvebin.points.pts1 = BSEDict["pts1"]
        _evolvebin.points.pts2 = BSEDict["pts2"]
        _evolvebin.points.pts3 = BSEDict["pts3"]
        _evolvebin.snvars.ecsn = BSEDict["ecsn"]
        _evolvebin.snvars.ecsn_mlow = BSEDict["ecsn_mlow"]
        _evolvebin.flags.aic = BSEDict["aic"]
        _evolvebin.ceflags.ussn = BSEDict["ussn"]
        _evolvebin.snvars.sigma = BSEDict["sigma"]
        _evolvebin.snvars.sigmadiv = BSEDict["sigmadiv"]
        _evolvebin.snvars.bhsigmafrac = BSEDict["bhsigmafrac"]
        _evolvebin.snvars.polar_kick_angle = BSEDict["polar_kick_angle"]
        _evolvebin.snvars.natal_kick_array = BSEDict["natal_kick_array"]
        _evolvebin.cevars.qcrit_array = BSEDict["qcrit_array"]
        _evolvebin.windvars.beta = BSEDict["beta"]
        _evolvebin.windvars.xi = BSEDict["xi"]
        _evolvebin.windvars.acc2 = BSEDict["acc2"]
        _evolvebin.windvars.epsnov = BSEDict["epsnov"]
        _evolvebin.windvars.eddfac = BSEDict["eddfac"]
        _evolvebin.windvars.gamma = BSEDict["gamma"]
        _evolvebin.flags.bdecayfac = BSEDict["bdecayfac"]
        _evolvebin.magvars.bconst = BSEDict["bconst"]
        _evolvebin.magvars.ck = BSEDict["ck"]
        _evolvebin.flags.windflag = BSEDict["windflag"]
        _evolvebin.flags.qcflag = BSEDict["qcflag"]
        _evolvebin.flags.eddlimflag = BSEDict["eddlimflag"]
        _evolvebin.tidalvars.fprimc_array = BSEDict["fprimc_array"]
        _evolvebin.rand1.idum1 = -1
        _evolvebin.flags.bhspinflag = BSEDict["bhspinflag"]
        _evolvebin.snvars.bhspinmag = BSEDict["bhspinmag"]
        _evolvebin.mixvars.rejuv_fac = BSEDict["rejuv_fac"]
        _evolvebin.flags.rejuvflag = BSEDict["rejuvflag"]
        _evolvebin.flags.htpmb = BSEDict["htpmb"]
        _evolvebin.flags.st_cr = BSEDict["ST_cr"]
        _evolvebin.flags.st_tide = BSEDict["ST_tide"]
        _evolvebin.snvars.rembar_massloss = BSEDict["rembar_massloss"]
        _evolvebin.metvars.zsun = BSEDict["zsun"]
        _evolvebin.snvars.kickflag = BSEDict["kickflag"]
        _evolvebin.cmcpass.using_cmc = 0

        # kstar, mass, orbital period (days), eccentricity, metaliccity, evolution time (millions of years)
        initial_stars = InitialBinaryTable.InitialBinaries(
            mass,
            np.ones_like(mass) * 0,
            np.ones_like(mass) * -1,
            np.ones_like(mass) * -1,
            np.ones_like(mass) * 0.1,
            self.set_kstar(mass),
            np.ones_like(mass) * 0,
            np.ones_like(mass) * metallicity,
        )
        initial_stars["dtp"] = initial_stars["tphysf"]

        initial_stars = initial_stars.assign(
            kick_info=[np.zeros((2, 17))] * len(initial_stars)
        )
        initial_conditions = initial_stars.to_dict("records")

        rad_1 = np.zeros(len(initial_stars))
        for idx, initial_condition in enumerate(initial_conditions):
            [bpp_index, bcm_index, _] = _evolvebin.evolv2(
                [initial_condition["kstar_1"], initial_condition["kstar_2"]],
                [initial_condition["mass_1"], initial_condition["mass_2"]],
                initial_condition["porb"],
                initial_condition["ecc"],
                initial_condition["metallicity"],
                initial_condition["tphysf"],
                initial_condition["dtp"],
                [initial_condition["mass0_1"], initial_condition["mass0_2"]],
                [initial_condition["rad_1"], initial_condition["rad_2"]],
                [initial_condition["lum_1"], initial_condition["lum_2"]],
                [initial_condition["massc_1"], initial_condition["massc_2"]],
                [initial_condition["radc_1"], initial_condition["radc_2"]],
                [initial_condition["menv_1"], initial_condition["menv_2"]],
                [initial_condition["renv_1"], initial_condition["renv_2"]],
                [initial_condition["omega_spin_1"], initial_condition["omega_spin_2"]],
                [initial_condition["B_1"], initial_condition["B_2"]],
                [initial_condition["bacc_1"], initial_condition["bacc_2"]],
                [initial_condition["tacc_1"], initial_condition["tacc_2"]],
                [initial_condition["epoch_1"], initial_condition["epoch_2"]],
                [initial_condition["tms_1"], initial_condition["tms_2"]],
                [initial_condition["bhspin_1"], initial_condition["bhspin_2"]],
                initial_condition["tphys"],
                np.zeros(20),
                np.zeros(20),
                initial_condition["kick_info"],
            )
            rad_1[idx] = _evolvebin.binary.bcm[0, 5]

        return rad_1
