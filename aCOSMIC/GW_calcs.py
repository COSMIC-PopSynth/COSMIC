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
from scipy.interpolate import interp1d 

__author__ = 'Katelyn Breivik <katie.breivik@gmail.com>'
__credits__ = 'Scott Coughlin <scott.coughlin@ligo.org>'
__all__ = 'GW_calcs'

G = 6.67384e-11
c = 2.99792458e8
parsec = 3.08567758e16
Rsun = 6.955e8
Msun = 1.9891e30
day = 86400.0
rsun_in_au = 215.0954
day_in_year = 365.242
sec_in_day = 86400.0
sec_in_hour = 3600.0
hrs_in_day = 24.0
sec_in_year = 3.15569e7
Tobs = 4*sec_in_year
geo_mass = G/c**2

def m_chirp(m1, m2):
    return (m1*m2)**(3./5.)/(m1+m2)**(1./5.)

def peak_gw_freq(m1, m2, ecc, porb):
    # UNITS ARE SI!!!
    # convert the orbital period into a separation using Kepler III
    sep_m = (G/(4*np.pi**2)*porb**2*(m1+m2))**(1./3.)
    
    f_gw_peak = ((G*(m1+m2))**0.5/np.pi) * (1+ecc)**1.1954/(sep_m*(1-ecc)**2)**1.5
    return f_gw_peak

def peters_gfac(ecc, n_harmonic):
    g_fac_squared = []    
    for n in range(1,n_harmonic):
        g_fac_squared.append((n**4 / 32.0)*( (ss.jv((n-2), (n*ecc)) - 2*ecc*ss.jv((n-1), (n*ecc)) +\
                             2.0/n*ss.jv((n), (n*ecc)) + 2*ecc*ss.jv((n+1), (n*ecc)) -\
                             ss.jv((n+2), (n*ecc)))**2 +\
                             (1-ecc**2)*(ss.jv((n-2), (n*ecc)) - 2*ss.jv((n), (n*ecc)) +\
                             ss.jv((n+2), (n*ecc)))**2 +\
                             (4/(3.0*n**2))*(ss.jv((n), (n*ecc)))**2
                             )/n**2.0)

    return np.array(g_fac_squared)
 
def LISA_calcs(m1, m2, porb, ecc, dist, n_harmonic):
    # UNITS ARE IN SI!!
    # LISA mission: 2.5 million km arms
    LISA_psd = lisa_sensitivity.lisa_psd()
    
    mChirp = m_chirp(m1, m2)
    h_0 = np.sqrt(32.0/5.0)*G/c**2*(mChirp)/(dist)*(G/c**3*np.pi*(1/(porb))*mChirp)**(2./3.)
    h_0_squared = 2*h_0**2
    #h_0_squared = (1.607e-21)**2 * ((mChirp)**(5.0/3.0) * (porb)**(-2.0/3.0) * dist**(-1.0))**2
    freq = []
    psd = []

    SNR = np.zeros(len(mChirp))
    gw_freq = np.zeros(len(mChirp))

    ind_ecc, = np.where(ecc > 0.1)
    ind_circ, = np.where(ecc <= 0.1)
    
    SNR[ind_circ] = (h_0_squared.iloc[ind_circ] * 1.0/4.0 * Tobs / LISA_psd(2 / porb.iloc[ind_circ]))
    gw_freq[ind_circ] = 2 / (porb.iloc[ind_circ])
    
    power = h_0_squared.iloc[ind_circ] * 1.0/4.0 * Tobs
    psd.extend(power)
    freq.extend(2.0 / (porb.iloc[ind_circ]))
    if len(ind_ecc) > 0:
        peters_g_factor = peters_gfac(ecc.iloc[ind_ecc], n_harmonic)
        GW_freq_array = np.array([np.arange(1,n_harmonic)/p for p in porb.iloc[ind_ecc]]).T
        GW_freq_flat = GW_freq_array.flatten()
        LISA_curve_eval_flat = LISA_psd(GW_freq_flat)
        LISA_curve_eval = LISA_curve_eval_flat.reshape((n_harmonic-1,len(ind_ecc)))

        h_0_squared_ecc = np.array([h*np.ones(n_harmonic-1) for h in h_0_squared.iloc[ind_ecc]]).T
        SNR_squared = h_0_squared_ecc * Tobs * peters_g_factor / LISA_curve_eval

        power = h_0_squared_ecc * Tobs * peters_g_factor
        power_flat = power.flatten()

        psd.extend(power_flat)
        freq.extend(GW_freq_flat)

        SNR[ind_ecc] = (SNR_squared.sum(axis=0))**0.5
        gw_freq[ind_ecc] = peak_gw_freq(m1.iloc[ind_ecc], m2.iloc[ind_ecc], ecc.iloc[ind_ecc], porb.iloc[ind_ecc])

    SNR_dat = pd.DataFrame(np.vstack([gw_freq, SNR, m1, m2, ecc, porb, dist]).T, columns=['gw_freq', 'SNR', 'm1', 'm2', 'ecc', 'porb', 'dist'])
    
    SNR_dat = SNR_dat.loc[SNR_dat.dist > 5.0*parsec]
    if len(SNR_dat.SNR > 1.0) > 1:
        SNR_dat = SNR_dat.loc[SNR_dat.SNR > 1.0]
    else:
        SNR_dat = pd.DataFrame(columns=['gw_freq', 'SNR', 'm1', 'm2', 'ecc', 'porb', 'dist']) 
        
    psd_dat = pd.DataFrame(np.vstack([freq, psd]).T, columns=['freq', 'PSD'])

    return SNR_dat, psd_dat.sort_values(by=['freq'])

def compute_foreground(psd_dat):
    binwidth = 40.0/Tobs
    freqBinsLISA = np.arange(1e-5,1e-1,binwidth)
    print 'number of bins in foreground', len(freqBinsLISA)
    binIndices = np.digitize(psd_dat.freq, freqBinsLISA)
    print 'bins digitized'
    psd_dat['digits'] = binIndices
    power_sum = psd_dat[['PSD', 'digits']].groupby('digits').sum()['PSD']
    power_tot = np.zeros(len(freqBinsLISA)+1)
    power_tot[power_sum.index[:len(freqBinsLISA)-1]] = power_sum
    foreground_dat = pd.DataFrame(np.vstack([freqBinsLISA, power_tot[1:]]).T,\
                                  columns=['freq', 'PSD'])
    return foreground_dat
