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
geo_mass = G/c**2

def m_chirp(m1, m2):
    """Computes the chirp mass in the units of mass supplied
   
    Parameters
    ----------
    m1 (float or array):
        primary mass
  
    m2 (float or array):    
        secondary mass

    Returns
    -------
    m_chirp (float or array):
        chirp mass in units of mass supplied
    """
    return (m1*m2)**(3./5.)/(m1+m2)**(1./5.)

def peak_gw_freq(m1, m2, ecc, porb):
    """Computes the peak gravitational-wave frequency for an
    eccentric binary system. Units are SI

    Paramters
    ---------
    m1 (float or array):
        primary mass [kg]
  
    m2 (float or array):    
        secondary mass [kg]

    ecc (float or array):
        eccentricity

    porb (float or array):
        orbital period [s]

    Returns
    -------
    f_gw_peak (floar or array):
        peak gravitational-wave frequency [Hz]
    """

    # convert the orbital period into a separation using Kepler III
    sep_m = (G/(4*np.pi**2)*porb**2*(m1+m2))**(1./3.)
    
    f_gw_peak = ((G*(m1+m2))**0.5/np.pi) * (1+ecc)**1.1954/(sep_m*(1-ecc)**2)**1.5
    return f_gw_peak

def peters_gfac(ecc, n_harmonic):
    """Computes the factor g_n/n**2 from Peters & Mathews 1963

    Parameters
    ----------
    ecc (float or array):
        eccentricity

    n_harmonic (int):
        number of frequency harmonics to include

    Returns
    -------
    g_fac_squared (array):
        array of g_n/n**2
    """
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
 
def LISA_SNR(m1, m2, porb, ecc, dist, n_harmonic, Tobs):
    """Computes the LISA signal-to-noise ratio with inputs in SI units
    with a LISA mission of 4 yr and 2.5 million km arms
    Note that this does not take into account the Galactic binary foreground
    and assumes no orbital frequency evolution due to GW emission

    Parameters
    ----------
    m1 (float or array):
        primary mass [kg]
  
    m2 (float or array):    
        secondary mass [kg]

    porb (float or array):
        orbital period [s]

    ecc (float or array):
        eccentricity

    dist (float or array):
        Solar-centric distance [m]

    n_harmonic (int):
        number of frequency harmonics to include

    Tobs (float):
            LISA observation time in seconds

    Returns
    -------
    SNR_dat (DataFrame):
        DataFrame with columns ['freq', 'SNR'] where
        freq = gravitational-wave frequency in Hz and 
        SNR = LISA signal-to-noise ratio 
    """

    # LISA mission: 2.5 million km arms
    LISA_psd = lisa_sensitivity.lisa_root_psd()
    
    mChirp = m_chirp(m1, m2)
    h_0 = 8*G/c**2*(mChirp)/(dist)*(G/c**3*2*np.pi*(1/(porb))*mChirp)**(2./3.)
    h_0_squared = h_0**2
  
    SNR = np.zeros(len(mChirp))
    gw_freq = np.zeros(len(mChirp))

    ind_ecc, = np.where(ecc > 0.1)
    ind_circ, = np.where(ecc <= 0.1)

    SNR[ind_circ] = (h_0_squared.iloc[ind_circ] * 1.0/4.0 * Tobs / LISA_psd(2 / porb.iloc[ind_circ])**2)
    gw_freq[ind_circ] = 2 / (porb.iloc[ind_circ])

    if len(ind_ecc) > 0:
        peters_g_factor = peters_gfac(ecc.iloc[ind_ecc], n_harmonic)
        GW_freq_array = np.array([np.arange(1,n_harmonic)/p for p in porb.iloc[ind_ecc]]).T
        GW_freq_flat = GW_freq_array.flatten()
        LISA_curve_eval_flat = LISA_psd(GW_freq_flat)**2
        LISA_curve_eval = LISA_curve_eval_flat.reshape((n_harmonic-1,len(ind_ecc)))

        h_0_squared_ecc = np.array([h*np.ones(n_harmonic-1) for h in h_0_squared.iloc[ind_ecc]]).T
        SNR_squared = h_0_squared_ecc * Tobs * peters_g_factor / LISA_curve_eval

        SNR[ind_ecc] = (SNR_squared.sum(axis=0))**0.5
        gw_freq[ind_ecc] = peak_gw_freq(m1.iloc[ind_ecc], m2.iloc[ind_ecc], ecc.iloc[ind_ecc], porb.iloc[ind_ecc])    
    SNR_dat = pd.DataFrame(np.vstack([gw_freq, SNR]).T, columns=['gw_freq', 'SNR'])
    if len(SNR_dat.SNR > 1.0) > 1:
        SNR_dat = SNR_dat.loc[SNR_dat.SNR > 1.0]
    else:
        SNR_dat = pd.DataFrame(columns=['gw_freq', 'SNR'])
    
    return SNR_dat


def LISA_PSD(m1, m2, porb, ecc, dist, n_harmonic, Tobs):
    """Computes the LISA power spectral density with inputs in SI units
    with a LISA mission of Tobs [sec] and 2.5 million km arms
    Note that this does not take into account the Galactic binary foreground
    and assumes no orbital frequency evolution due to GW emission

    Parameters
    ----------
    m1 (float or array):
        primary mass [kg]
  
    m2 (float or array):    
        secondary mass [kg]

    porb (float or array):
        orbital period [s]

    ecc (float or array):
        eccentricity

    dist (float or array):
        Solar-centric distance [m]

    n_harmonic (int):
        number of frequency harmonics to include

    Tobs (float):
            LISA observation time in seconds

    Returns
    -------
    PSD_dat (DataFrame):
        DataFrame with columns ['freq', 'PSD'] where
        freq = gravitational wave frequency in Hz
        SNR = LISA signal-to-noise ratio 
    """

    # LISA mission: 2.5 million km arms
    LISA_psd = lisa_sensitivity.lisa_root_psd()
    
    mChirp = m_chirp(m1, m2)
    h_0 = 8*G/c**2*(mChirp)/(dist)*(G/c**3*2*np.pi*(1/(porb))*mChirp)**(2./3.)
    h_0_squared = h_0**2

    freq = []
    psd = []

    ind_ecc, = np.where(ecc > 0.1)
    ind_circ, = np.where(ecc <= 0.1)

    power = h_0_squared.iloc[ind_circ] * 1.0/4.0 * Tobs
    
    psd.extend(power)
    freq.extend(2.0 / (porb.iloc[ind_circ]))

    if len(ind_ecc) > 0:
        peters_g_factor = peters_gfac(ecc.iloc[ind_ecc], n_harmonic)
        power = h_0_squared_ecc * Tobs * peters_g_factor
        power_flat = power.flatten()

        psd.extend(power_flat)
        freq.extend(GW_freq_flat)
    PSD_dat = pd.DataFrame(np.vstack([freq, psd]).T, columns=['freq', 'PSD'])
    
    return PSD_dat 

def compute_foreground(psd_dat, Tobs=4*sec_in_year):
    """Computes the gravitational-wave foreground by binning the PSDs 
    in psd_dat according to the LISA frequency resolution where 1
    frequency bin has a binwidth of 1/(Tobs [s])

    Parameters
    ----------
    psd_dat (DataFrame):
        DataFrame with columns ['freq', 'PSD'] where 
        freq = gravitational-wave frequency in Hz
        PSD = LISA power spectral density

    Tobs (float):
        LISA observation time in seconds; Default=4 yr

    Returns
    -------
    foreground_dat (DataFrame):
        DataFrame with columns ['freq', 'PSD'] where 
        freq = gravitational-wave frequency of LISA frequency bins in Hz
        PSD = LISA power spectral density at each LISA frequency
    """

    binwidth = 1.0/Tobs
    freqBinsLISA = np.arange(5e-4,1e-2,binwidth)
    binIndices = np.digitize(psd_dat.freq, freqBinsLISA)
    psd_dat['digits'] = binIndices
    power_sum = psd_dat[['PSD', 'digits']].groupby('digits').sum()['PSD']
    power_tot = np.zeros(len(freqBinsLISA)+1)
    power_tot[power_sum.index[:len(freqBinsLISA)-1]] = power_sum
    foreground_dat = pd.DataFrame(np.vstack([freqBinsLISA, power_tot[1:]]).T,\
                                  columns=['freq', 'PSD'])
    return foreground_dat
