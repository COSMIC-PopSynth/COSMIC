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
    # LISA mission: 2.5 million km arms
    LISA_root_psd = lisa_sensitivity.lisa_root_psd(np.logspace(-7, 1, 1000), 2.5e9)
    h_0_squared = 1.0e-42 * ((mChirpList)**(5.0/3.0) * (porbList)**(-2.0/3.0) * distList**(-1.0))**2
    
    ind_ecc, = np.where(eccList > 0.1)
    ind_circ, = np.where(eccList <= 0.1)
   
    SNR = np.zeros(len(mChirpList))
    
    SNR[ind_circ] = (h_0_squared[ind_circ] * Tobs / LISA_root_psd(2 / (porbList[ind_circ] * sec_in_hour))**2)**0.5
    
    for ecc, porb, ind in zip(eccList[ind_ecc], porbList[ind_ecc], ind_ecc):
        SNR_squared = np.zeros(n_harmonic)
        for n in range(1, n_harmonic):
            SNR_squared[n] = h_0_squared[ind] * peters_gfac(ecc, n) *\
                             Tobs / LISA_root_psd(n / (porb * sec_in_hour))**2
        SNR[ind] = np.sum(SNR_squared)**0.5
    
    return SNR

def power_spectral_density(mChirpList, porbList, eccList, distList, n_harmonic):
    psd = []
    freq = []
    
    # chirp mass in Msun, porb in hr, dist in kpc
    h_0_squared = 1.0e-42 * ((mChirpList)**(5.0/3.0) * (porbList)**(-2.0/3.0) * distList**(-1.0))**2

    ind_ecc, = np.where(eccList > 0.1)
    ind_circ, = np.where(eccList <= 0.1)

    psd.extend(h_0_squared[ind_circ] * Tobs)
    freq.extend(2.0 / (porbList[ind_circ] * sec_in_hour))

    for ecc, porb, ind in zip(eccList[ind_ecc], porbList[ind_ecc], ind_ecc):
        psd.extend(h_0_squared[ind] * peters_gfac(ecc, np.arange(1,n_harmonic)) * Tobs)
        freq.extend(np.arange(1,n_harmonic) / (porb * sec_in_hour))
    psd_dat = pd.DataFrame(np.vstack([freq, psd]).T, columns=['freq', 'psd'])
    return psd_dat.sort_values(by=['freq'])
    
def moving_average(interval, window_size):
    window = np.ones(int(window_size))/float(window_size)
    return np.convolve(interval, window, 'same')

def compute_foreground(psd_dat):
    binwidth = 4.0/Tobs
    freqBinsLISA = np.arange(1e-6,1e-1,binwidth)
    print 'number of bins in foreground', len(freqBinsLISA)
    binIndices = np.digitize(psd_dat.freq, freqBinsLISA)
    print 'bins digitized'
    psd_dat['digits'] = binIndices
    power_sum = psd_dat[['psd', 'digits']].groupby('digits').sum()['psd']
    power_tot = np.zeros(len(freqBinsLISA))
    power_tot[power_sum.index[:len(freqBinsLISA)]] = power_sum
    foreground_dat = pd.DataFrame(np.vstack([freqBinsLISA, power_tot]).T,\
                                  columns=['freq', 'psd'])
    return foreground_dat
