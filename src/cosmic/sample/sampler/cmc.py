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

"""`cmc`
"""

import numpy as np

from .sampler import register_sampler
from .independent import Sample
from .. import InitialCMCTable, InitialBinaryTable
from ..cmc import elson, king
from ... import utils

__author__ = "Scott Coughlin <scott.coughlin@ligo.org>"
__credits__ = [
   "Carl Rodriguez <carllouisrodriguez@gmail.com>", 
   "Newlin Weatherford <newlinweatherford2017@u.northwestern.edu>"
]
__all__ = ["get_cmc_sampler", "CMCSample"]


def get_cmc_sampler(
    cluster_profile, primary_model, ecc_model, porb_model, binfrac_model, met, size, **kwargs
):
    """Generates an initial cluster sample according to user specified models

    Parameters
    ----------
    cluster_profile : `str`
        Model to use for the cluster profile (i.e. sampling of the placement of objects in the cluster and their velocity within the cluster)
        Options include:
        'plummer' : Standard Plummer sphere.
            Additional parameters:
            'r_max' : `float`
                the maximum radius (in virial radii) to sample the clsuter
        'elson' : EFF 1987 profile.  Generalization of Plummer that better fits young massive clusters
            Additional parameters:
            'gamma' : `float`
                steepness paramter for Elson profile; note that gamma=4 is same is Plummer
            'r_max' : `float`
                the maximum radius (in virial radii) to sample the clsuter
        'king' : King profile 
            'w_0' : `float`
                King concentration parameter
            'r_max' : `float`
                the maximum radius (in virial radii) to sample the clsuter

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

    binfrac_model : `str or float`
        Model for binary fraction; choices include: vanHaaften, offner22, or a fraction where 1.0 is 100% binaries

    binfrac_model_msort : `str or float`
        Same as binfrac_model for M>msort

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

    met : `float`
        Sets the metallicity of the binary population where solar metallicity is zsun 

    size : `int`
        Size of the population to sample

    zsun : `float`
        optional kwarg for setting effective radii, default is 0.02

    Optional Parameters
    -------------------
    virial_radius : `float`
        the initial virial radius of the cluster, in parsecs
        Default -- 1 pc

    tidal_radius : `float`
        the initial tidal radius of the cluster, in units of the virial_radius
        Default -- 1e6 rvir

    central_bh : `float`
        Put a central massive black hole in the cluster
        Default -- 0 MSUN
        
    scale_with_central_bh : `bool`
        If True, then the potential from the central_bh is included when scaling the radii.
        Default -- False

    seed : `float`
        seed to the random number generator, for reproducability

    Returns
    -------
    Singles: `pandas.DataFrame`
        DataFrame of Single objects in the format of the InitialCMCTable

    Binaries: `pandas.DataFrame`
        DataFrame of Single objects in the format of the InitialCMCTable
    """
    initconditions = CMCSample()

    # if RNG seed is provided, then use it globally
    rng_seed = kwargs.pop("seed", 0)
    if rng_seed != 0:
        np.random.seed(rng_seed)

    # get radii, radial and transverse velocities
    r, vr, vt = initconditions.set_r_vr_vt(cluster_profile, N=size, **kwargs)

    # track the mass in singles and the mass in binaries
    mass_singles = 0.0
    mass_binaries = 0.0

    # track the total number of stars sampled
    n_singles = 0
    n_binaries = 0

    mass1, total_mass1 = initconditions.sample_primary(
        primary_model, size=size, **kwargs)
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

    # set singles id
    single_ids = np.arange(mass1.size)
    binary_secondary_object_id = np.arange(mass1.size, mass1.size + mass2_binaries.size)

    # set binind and correct masses of binaries
    binind = np.zeros(mass1.size)
    binind[binary_index] = np.arange(len(binary_index)) + 1
    mass1[binary_index] += mass2_binaries

    # track the total number sampled
    n_singles += len(mass_single)
    n_binaries += len(mass1_binaries)

    # Obtain radii (technically this is done for the binaries in the independent sampler 
    # if set_radii_with_BSE is true, but that's not a huge amount of overhead)
    zsun = kwargs.pop("zsun", 0.02)

    Reff = initconditions.set_reff(mass1, metallicity=met, zsun=zsun)
    Reff1 = Reff[binary_index]
    Reff2 = initconditions.set_reff(mass2_binaries, metallicity=met, zsun=zsun)

    # select out the primaries and secondaries that will produce the final kstars
    porb_max = initconditions.calc_porb_max(mass1, vr, vt, binary_index, mass1_binaries, mass2_binaries, **kwargs)

    porb,aRL_over_a = initconditions.sample_porb(
        mass1_binaries, mass2_binaries, Reff1, Reff2, porb_model=porb_model, porb_max=porb_max, size=mass1_binaries.size
    )
    ecc = initconditions.sample_ecc(aRL_over_a, ecc_model, size=mass1_binaries.size)

    sep = utils.a_from_p(porb, mass1_binaries, mass2_binaries)
    kstar1 = initconditions.set_kstar(mass1_binaries)
    kstar2 = initconditions.set_kstar(mass2_binaries)

    singles_table = InitialCMCTable.InitialCMCSingles(
        single_ids +
        1, initconditions.set_kstar(mass1), mass1, Reff, r, vr, vt, binind
    )

    binaries_table = InitialCMCTable.InitialCMCBinaries(
        np.arange(mass1_binaries.size) + 1,
        single_ids[binary_index] + 1,
        kstar1,
        mass1_binaries,
        Reff1,
        binary_secondary_object_id + 1,
        kstar2,
        mass2_binaries,
        Reff2,
        sep,
        ecc,
    )

    singles_table.metallicity = met
    binaries_table.metallicity = met
    singles_table.virial_radius = kwargs.get("virial_radius",1) 
    singles_table.tidal_radius = kwargs.get("tidal_radius",1e6) 
    singles_table.central_bh = kwargs.get("central_bh",0)
    singles_table.scale_with_central_bh = kwargs.get("scale_with_central_bh",False)
    singles_table.mass_of_cluster = np.sum(singles_table["m"]) + singles_table.central_bh

    return singles_table, binaries_table


register_sampler(
    "cmc",
    InitialCMCTable,
    get_cmc_sampler,
    usage="primary_model ecc_model porb_model binfrac_model met size",
)

def get_cmc_point_mass_sampler(
    cluster_profile, size, **kwargs
):
    """Generates an CMC cluster model according to user-specified model.
    Note here that masses will all be unity (with total cluster normalized accordingly)

    Parameters
    ----------
    cluster_profile : `str`
        Model to use for the cluster profile (i.e. sampling of the placement of objects in the cluster and their velocity within the cluster)
        Options include:
        'plummer' : Standard Plummer sphere.
            Additional parameters:
            'r_max' : `float`
                the maximum radius (in virial radii) to sample the clsuter
        'elson' : EFF 1987 profile.  Generalization of Plummer that better fits young massive clusters
            Additional parameters:
            'gamma' : `float`
                steepness paramter for Elson profile; note that gamma=4 is same is Plummer
            'r_max' : `float`
                the maximum radius (in virial radii) to sample the clsuter
        'king' : King profile 
            'w_0' : `float`
                King concentration parameter
            'r_max' : `float`
                the maximum radius (in virial radii) to sample the clsuter

    size : `int`
        Size of the population to sample

    Optional Parameters
    -------------------
    virial_radius : `float`
        the initial virial radius of the cluster, in parsecs
        Default -- 1 pc

    tidal_radius : `float`
        the initial tidal radius of the cluster, in units of the virial_radius
        Default -- 1e6 rvir

    central_bh : `float`
        Put a central massive black hole in the cluster.  Note here units are 
        in code units where every particle has mass 1 (before normalization), e.g. if you
        want a cluster with 10000 stars and a BH that is 10% of the total mass, central_bh=1000
        Default -- 0
        
    scale_with_central_bh : `bool`
        If True, then the potential from the central_bh is included when scaling the radii.
        Default -- False

    seed : `float`
        seed to the random number generator, for reproducability

    Returns
    -------
    Singles: `pandas.DataFrame`
        DataFrame of Single objects in the format of the InitialCMCTable
    Binaries: `pandas.DataFrame`
        DataFrame of Single objects in the format of the InitialCMCTable

    """
    initconditions = CMCSample()

    # if RNG seed is provided, then use it globally
    rng_seed = kwargs.pop("seed", 0)
    if rng_seed != 0:
        np.random.seed(rng_seed)

    # get radii, radial and transverse velocities
    r, vr, vt = initconditions.set_r_vr_vt(cluster_profile, N=size, **kwargs)

    mass1 = np.ones(size)/size
    Reff = np.zeros(size)
    binind = np.zeros(size)
    single_ids = np.arange(size)

    singles_table = InitialCMCTable.InitialCMCSingles(
        single_ids + 1, initconditions.set_kstar(mass1), mass1, Reff, r, vr, vt, binind
    )

    # We assume no binaries for the point mass models
    binaries_table = InitialCMCTable.InitialCMCBinaries(1,1,0,0,0,1,0,0,0,0,0)

    singles_table.metallicity = 0.02
    binaries_table.metallicity = 0.02
    singles_table.virial_radius = kwargs.get("virial_radius",1)
    singles_table.tidal_radius = kwargs.get("tidal_radius",1e13)
    singles_table.central_bh = kwargs.get("central_bh",0)
    singles_table.scale_with_central_bh = kwargs.get("scale_with_central_bh",False)
    singles_table.mass_of_cluster = np.sum(singles_table["m"])*size + singles_table.central_bh

    # Already scaled from the IC generators (unless we've added a central BH)
    if singles_table.central_bh != 0:
        singles_table.scaled_to_nbody_units = False
    else:
        singles_table.scaled_to_nbody_units = True

    # Already scaled from the IC generators (unless we've added a central BH)
    if singles_table.central_bh != 0:
        singles_table.scaled_to_nbody_units = False 
    else:
        singles_table.scaled_to_nbody_units = True
   
    return singles_table, binaries_table

register_sampler(
    "cmc_point_mass",
    InitialCMCTable,
    get_cmc_point_mass_sampler,
    usage="size",
)

class CMCSample(Sample):
    def set_r_vr_vt(self,cluster_profile, **kwargs):
        if cluster_profile == "elson":
            elson_kwargs = {
                k: v for k, v in kwargs.items() if k in ["gamma", "r_max", "N"]
            }
            r, vr, vt = elson.draw_r_vr_vt(**elson_kwargs)
        elif cluster_profile == "plummer":
            plummer_kwargs = {k: v for k, v in kwargs.items() if k in ["r_max", "N"]}
            r, vr, vt = elson.draw_r_vr_vt(gamma=4, **plummer_kwargs)
        elif cluster_profile == "king":
            king_kwargs = {k: v for k, v in kwargs.items() if k in ["w_0", "N"]}
            r, vr, vt = king.draw_r_vr_vt(**king_kwargs)
        else:
            raise ValueError("Cluster profile not defined! Please specify one of [elson, plummer, king]")

        return r, vr, vt

    def calc_porb_max(self, mass, vr, vt, binary_index, mass1_binary, mass2_binary, **kwargs): 

        ## First, compute a rolling average with window length of AVEKERNEL
        ## NOTE: this is in cluster code units (i.e. M=G=1), same as vr and vt
        AVEKERNEL = 20

        ## Then compute the average mass and avergae of m*v^2
        ## First, to do this properly, we need to reflect the boundary points (so that the average 
        ## at the boundary doesn't go to zero artifically))
        m = np.concatenate([mass[AVEKERNEL-1::-1],mass,mass[:-AVEKERNEL-1:-1]])
        mv2 = mass*(vr*vr + vt*vt) 
        mv2 = np.concatenate([mv2[AVEKERNEL-1::-1],mv2,mv2[:-AVEKERNEL-1:-1]])

        ## Then compute the averages by convolving with a uniform array (note we trim the reflected points off here)
        m_average = np.convolve(m,np.ones(AVEKERNEL),mode='same')[AVEKERNEL:-AVEKERNEL]/AVEKERNEL
        mv2_average = np.convolve(mv2,np.ones(AVEKERNEL),mode='same')[AVEKERNEL:-AVEKERNEL]/AVEKERNEL

        ## Finally return the mass-weighted average 3D velocity dispersion (in cluster code units)
        sigma =  np.sqrt(mv2_average/m_average)

        ## Now compute the orbital velocity corresponding to the hard/soft boundary; we only want it for the binaries
        v_orb = 0.7*1.30294*sigma[binary_index]   #sigma * 4/sqrt(3/pi); 0.7 is a factor we use in CMC
        #v_orb = 0.7*1.30294*np.mean(sigma[:20])   #old CMC way of using just core velocity dispersion; TODO: possibly make flag option?


        ## Maximum semi-major axis just comes from Kepler's 3rd
        ## Note, to keep in code units, we need to divide binary masses by total cluster mass
        amax = (mass1_binary+mass2_binary) / v_orb**2 / np.sum(mass) 

        ## Convert from code units (virial radii) to RSUN 
        virial_radius = kwargs.get("virial_radius",1) ## get the virial radius of the cluster (uses 1pc if not given)
        RSUN_PER_PARSEC = 4.435e+7
        amax *= RSUN_PER_PARSEC * virial_radius 

        ## Finally go from sep to porb
        porb_max = utils.p_from_a(amax, mass1_binary, mass2_binary)

        return porb_max ## returns orbital period IN DAYS


