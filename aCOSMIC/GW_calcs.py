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

'''GW_calcs
'''

import numpy as np
import math
import random
import scipy.special as ss
import scipy.stats as stats
import lisa_sensitivity
import pandas as pd

__author__ = 'Katelyn Breivik <katie.breivik@gmail.com>'
__credits__ = 'Scott Coughlin <scott.coughlin@ligo.org>'
__all__ = 'GW_calcs'

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
Tobs = 4*sec_in_year
geo_mass = G/c**2

def m_chirp(m1, m2):
    return (m1*m2)**(3./5.)/(m1+m2)**(1./5.)

def peters_gfac(ecc, n):
    # Note this is: g(n,e)/n**2 that is being calculated
    g_fac_squared = (n**4 / 32.0)*( (ss.jv((n-2), (n*ecc)) - 2*ecc*ss.jv((n-1), (n*ecc)) +\
                        2.0/n*ss.jv((n), (n*ecc)) + 2*ecc*ss.jv((n+1), (n*ecc)) -\
                        ss.jv((n+2), (n*ecc)))**2 +\
                       (1-ecc**2)*(ss.jv((n-2), (n*ecc)) - 2*ss.jv((n), (n*ecc)) +\
                       ss.jv((n+2), (n*ecc)))**2 +\
                       (4/(3.0*n**2))*(ss.jv((n), (n*ecc)))**2
                       )/n**2.0

    return g_fac_squared

def LISA_SNR(mChirpList, porbList, eccList, distList, n_harmonic):
    # chirp mass in Msun, porb in hr, dist in kpc
    LISA_root_psd = lisa_sensitivity.lisa_sensitivity()
    h_0_squared = 1.0e-42 * ((mChirpList)**(5.0/3.0) * (porbList)**(-2.0/3.0) * distList**(-1.0))**2
    
    if eccList.max() < 0.01:
        SNR_squared = h_0_squared * Tobs / LISA_root_psd(n / (porb * sec_in_hour))**2
        SNR_array = SNR_squared**0.5
    
    else:
        SNR_array = np.zeros(len(mChirpList))
        ii = 0
        for mChirp, porb, ecc, dist in zip(mChirpList, porbList, eccList, distList):
            SNR_squared = np.zeros(n_harmonic)
            for n in range(1, n_harmonic):
                SNR_squared[n] = h_0_squared[ii] * peters_gfac(ecc, n) *\
                                 Tobs / LISA_root_psd(n / (porb * sec_in_hour))**2
                SNR_array[ii] = np.sum(SNR_squared)**0.5
            ii += 1
    return SNR_array

def power_spectral_density(mChirpList, porbList, eccList, distList, n_harmonic):
    # chirp mass in Msun, porb in hr, dist in kpc
    psd = []
    freq = []

    # chirp mass in Msun, porb in hr, dist in kpc
    LISA_root_psd = lisa_sensitivity.lisa_sensitivity()
    h_0_squared = 1.0e-42 * ((mChirpList)**(5.0/3.0) * (porbList)**(-2.0/3.0) * distList**(-1.0))**2

    if eccList.max() < 0.01:
        psd = h_0_squared * Tobs / LISA_root_psd(2 / (porbList * sec_in_hour))**2
        freq = 2.0 / (porbList * sec_in_hour)
                                        
    else:
        SNR_array = np.zeros(len(mChirpList))
        for mChirp, porb, ecc, dist in zip(mChirpList, porbList, eccList, distList):
            for n in range(1, n_harmonic):
                psd.append(1.0e-42 * peters_gfac(ecc, n) * Tobs *\
                           ((mChirp)**(5.0/3.0) * (porb)**(-2.0/3.0) * dist**(-1.0))**2)
                freq.append(n / (porb * sec_in_hour))
    psd_dat = pd.DataFrame(np.vstack([freq, psd]).T, columns=['freq', 'psd'])
    return psd_dat
    
def moving_average(interval, window_size):
    window = np.ones(int(window_size))/float(window_size)
    return np.convolve(interval, window, 'same')

def compute_foreground(forb, power):
    nBinsLISA = 5000
    freqBinsLISA = np.logspace(-6,-2.5,nBinsLISA)
    print 'number of bins in foreground', nBinsLISA
    binIndices = np.digitize(forb,freqBinsLISA)
    print 'bins digitized'
    power_tot = [power[binIndices == ii].sum() for ii in range(1, len(freqBinsLISA))]
    
    #binValue = []
    #for ii in range(len(freqBinsLISA)):
    #    binIndex, =np.where(binIndices==ii)
    #    if len(binIndex) > 0:
    #        binValue.append(sum(power[binIndex]))
    #    else:
    #        binValue.append(0.0)

    print np.shape(power_tot), np.shape(freqBinsLISA)
    foreground_dat = pd.DataFrame(np.vstack([freqBinsLISA[1:], power_tot]).T,\
                                  columns=['freq', 'psd'])
    return foreground_dat
