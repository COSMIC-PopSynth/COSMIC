# -*- coding: utf-8 -*-
# Copyright (C) Scott Coughlin (2017)
#
# This file is part of aCOSMIC.
#
# aCOSMIC is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# aCOSMIC is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with aCOSMIC.  If not, see <http://www.gnu.org/licenses/>.

'''MC_samp
'''

import numpy as np
import math
import random
import scipy.special as ss
import scipy.stats as stats

__author__ = 'Katelyn Breivik <katie.breivik@gmail.com>'
__credits__ = 'Scott Coughlin <scott.coughlin@ligo.org>'
__all__ = 'MC_samp'


G = 6.67384*math.pow(10, -11.0)
c = 2.99792458*math.pow(10, 8.0)
parsec = 3.08567758*math.pow(10, 16)
Rsun = 6.955*math.pow(10, 8)
Msun = 1.9891*math.pow(10,30)
day = 86400.0
rsun_in_au = 215.0954
day_in_year = 365.242
sec_in_day = 86400.0
sec_in_hour = 3600.0
hrs_in_day = 24.0
sec_in_year = 3.15569*10**7.0
Tobs = 3.15569*10**7.0
geo_mass = G/c**2

# solar coordinates in the galaxy: in parsecs from 
# (Chaper 8 of Galactic Structure and stellar Pops book) Yoshii (2013)
############################################################################
x_sun = 8000.0
y_sun = 0.0
z_sun = 25

def mass_weighted_number(dat, total_sampled_mass, component_mass): 
    '''
    Compute the total number of systems in the synthetic catalog
    based on the total sampled mass of the simulated system and
    the total mass of a given galactic component
    '''
    print 'N binaries Fixed: ',dat.mass_1.size
    print 'Gx component mass: ',component_mass
    print 'Total sampled mass: ',total_sampled_mass 
      
    nSystems = int(dat.size*component_mass/total_sampled_mass)
    print 'N binaries in Gx: ',nSystems
    return nSystems

def select_component_mass(gx_component):
    # SET GALACTIC COMPONENT MASS ACCORDING TO MCMILLAN 2011
    # FOLLOWING BEST FIT
    ###########################################################################
    if gx_component == 'ThinDisk':
        gx_component_mass = 4.32e10
    elif gx_component == 'Bulge':
        gx_component_mass = 8.9e9
    elif gx_component == 'ThickDisk':
        gx_component_mass = 1.44e10
 
    return gx_component_mass

def sample_sech_squared(size, scale_height=0.3):
    # SAMPLE FROM SECH^2 DIST W/ PROVIDED SCALE HEIGHT FROM -INF to INF
    ###########################################################################
    rand_nums = np.random.uniform(0, 1, size)
    distributed_nums = np.arctanh( 2*rand_nums ) * scale_height

    return distributed_nums

def sample_exponential_radial(size, scale_height):
    # SAMPLE FROM EXP DIST W/ PROVIDED SCALE HEIGHT FROM 0 to INF
    ###########################################################################
    rand_nums = np.random.uniform(0, 1, size)
    distributed_nums = scale_height * np.log(1.0 - rand_nums)

    return distributed_nums

def sample_exponential_vertical(size, scale_height):
    # SAMPLE FROM EXP DIST W/ PROVIDED SCALE HEIGHT FROM 0 to INF
    ###########################################################################
    rand_nums = np.random.uniform(0, 1, size)
    distributed_nums = scale_height * np.log(1.0 - rand_nums)

    # CHOOSE A POSITION ABOVE AND BELOW THE DISK 
    ###########################################################################
    pos_neg_choose = np.random.uniform(0, 1, size)
    negInd, = np.where(pos_neg_choose < 0.5)
    distributed_nums[negInd] = distributed_nums[negInd] * (-1.0)

    return distributed_nums

def sample_exponential_square(size, scale_height):
    # SAMPLE FROM EXP(-(X/A)^2) DIST W/ PROVIDED SCALE HEIGHT FROM 0 to INF
    ###########################################################################
    rand_nums = np.random.uniform(0, 1, size)
    distributed_nums = scale_height * ss.erfinv(rand_nums)

    return distributed_nums

def galactic_position_sample(gx_component, size, model):
    # SAMPLE FROM VERY SIMPLE DISTRIBUTION FUNCTIONS 
    ###########################################################################
    if gx_component == 'ThinDisk':
        # COMPUTE THE POSITION AND ORIENTATION OF EACH BINARY
        #######################################################################
        # Assign the radial and vertical positions:
        if model == 'double_exp':
            r = sample_exponential_radial(size, 2.5)
            z = sample_exponential_vertical(size, 0.3)

        if model == 'sech_squared':
            r = sample_exponential_radial(size, 2.5)
            z = sample_sech_squared(size, 0.3)

        if model == 'McMillan':
            r = sample_exponential_radial(size, 2.9)
            z = sample_exponential_vertical(size, 0.3)

        # Assign the azimuthal positions:
        phi = np.random.uniform(0, 2*np.pi, size)
        
        # convert to cartesian and parsecs:
        xGX = r * np.cos(phi) * 1000.0
        yGX = r * np.sin(phi) * 1000.0
        zGX = z * 1000.0
    
    elif gx_component == 'Bulge':

        # COMPUTE THE POSITION AND ORIENTATION OF EACH BINARY
        #######################################################################
        # Assign the radial positions:
        if model == 'exp_squared':
            r = sample_exponential_square(size, 0.5)
            # Assign the polar positions (uniform in cos(theta):
            theta = np.pi - np.arccos(np.random.uniform(-1, 1, size))
        
            # Assign the azimuthal positions:
            phi = np.random.uniform(0, 2*np.pi , size)

            # convert to cartesian and parsecs:
            xGX = r * np.cos(phi) * np.sin(theta) * 1000.0
            yGX = r * np.sin(phi) * np.sin(theta) * 1000.0
            zGX = r * np.cos(theta) * 1000.0
        elif model == 'McMillan':
            r_save = []
            z_save = []
            # sample double exp func and then rejection sample
            while len(r_save) < size:
                rcut = 2.1
                q = 0.5
                r0 = 0.075
                alpha = -1.8
                r = np.random.uniform(0,1,size)
                z = np.random.uniform(0,1,size)
                prob = np.random.uniform(0,1,size)
                sample_func = np.exp(-(r**2 + (z/q)**2)/rcut**2)
                actual_func = (1+np.sqrt((r**2 + (z/q)**2))/r0)**(alpha)*sample_func
                indSave, = np.where(prob < actual_func)
                for ii in indSave:
                    r_save.append(r[ii])
                    z_save.append(z[ii])
                print 'length of sample: ',len(r_save)
            r = np.array(r_save[:size])
            z = np.array(z_save[:size])

            # Assign the azimuthal positions:
            phi = np.random.uniform(0, 2*np.pi , size)
            
            # convert to cartesian and parsecs:
            xGX = r * np.cos(phi) * 1000.0
            yGX = r * np.sin(phi) * 1000.0
            zGX = z * 1000.0

    elif gx_component == 'ThickDisk':
        # COMPUTE THE POSITION AND ORIENTATION OF EACH BINARY
        ###############################################################################
        # Assign the radial and vertical positions:
        if model == 'double_exp':
            r = sample_exponential_radial(size, 2.5)
            z = sample_exponential_vertical(size, 1.0)
        
        if model == 'McMillan':
            r = sample_exponential_radial(size, 3.31)
            z = sample_exponential_vertical(size, 0.9)

        # Assign the azimuthal position of the star
        phi = np.random.uniform(0, 2*np.pi, size)

        # convert to cartesian and parsecs
        xGX = r*np.cos(phi)*1000.0
        yGX = r*np.sin(phi)*1000.0
        zGX = z*1000.0


    # compute the distance to Earth/LISA/us in kiloparsecs
    dist = ((xGX- x_sun)**2 + (yGX - y_sun)**2 + (zGX - z_sun)**2)**(1/2.0)
    dist_kpc = dist/1000.0

    # assign an inclination, argument of periapsis, and longitude of ascending node
    inc = np.pi - np.arccos(np.random.uniform(-1, 0, size))
    OMEGA = np.random.uniform(0, 2*np.pi, size)
    omega = np.random.uniform(0, 2*np.pi, size)

         
    return np.vstack([xGX, yGX, zGX, dist_kpc, inc, OMEGA, omega])


def sample(dat, total_sampled_mass, gx_component, model='McMillan'):
    '''
    Monte-Carlo samples a galactic realization
    '''
    component_mass = select_component_mass(gx_component)
    nSystems = mass_weighted_number(dat, total_sampled_mass, component_mass)
    dat_trans = dat_transform(dat)
    
    # FIRST SAMPLE THE BINARY PARAMETERS FROM THE FIXED POPULATION
    # BY GENERATING A KDE FIT
    #####################################################################
    dat_kernel = stats.gaussian_kde(dat_trans)
    binary_sample_dat_trans = dat_kernel.resample(nSystems)
    binary_sample_dat = dat_un_transform(binary_sample_dat_trans, dat)
    
    # NEXT SAMPLE THE BINARY POSITIONS AND ORIENTATIONS
    #####################################################################
    binary_sample_positions = galactic_position_sample(gx_component, size = nSystems, model=model) 

    full_sample = np.concatenate([binary_sample_dat, binary_sample_positions])
    return full_sample
