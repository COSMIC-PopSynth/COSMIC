# -*- coding: utf-8 -*-
# Copyright (C) Scott Coughlin (2017)
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
__all__ = ['m_chirp', 'peak_gw_freq', 'peters_gfac', 'LISA_SNR', 'LISA_PSD', 'compute_foreground']

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
x_sun = 8.5
y_sun = 0.0
z_sun = 0.02

def m_chirp(m1, m2):
    """Computes the chirp mass in the units of mass supplied
   
    Parameters
    ----------
    m1 : float or array
        primary mass
    m2 : float or array    
        secondary mass

    Returns
    -------
    mchirp : float or array
        chirp mass in units of mass supplied
    """
    return (m1*m2)**(3./5.)/(m1+m2)**(1./5.)

def peak_gw_freq(m1, m2, ecc, porb):
    """Computes the peak gravitational-wave frequency for an
    eccentric binary system. Units are SI

    Parameters
    ----------
    m1 : float or array
        primary mass [kg]
    m2 : float or array    
        secondary mass [kg]
    ecc : float or array
        eccentricity
    porb : float or array
        orbital period [s]

    Returns
    -------
    f_gw_peak : float or array
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
    ecc : float or array
        eccentricity
    n_harmonic : int
        number of frequency harmonics to include

    Returns
    -------
    g_fac_squared : array
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
 
def LISA_combo(foreground_dat):
    """Computes the combination of the population foreground and LISA PSD
    
    Parameters
    ----------
    foreground_dat : DataFrame
        gw_freq [Hz], PSD [Hz^-1]

    Returns
    -------
    LISA_curve_combo : interpolation
        gw_freq [Hz], PSD [Hz^-1]
    """
    # LISA mission: 2.5 million km arms
    LISA_psd = lisa_sensitivity.lisa_root_psd()
    LISA_curve = LISA_psd(np.array(foreground_dat.freq))**2
    
    LISA_curve_foreground = LISA_curve + np.array(foreground_dat.PSD)
    LISA_combo_interp = interp1d(np.array(foreground_dat.freq), LISA_curve_foreground)
    
    return LISA_combo_interp


def LISA_SNR(m1, m2, porb, ecc, xGx, yGx, zGx, n_harmonic, Tobs, foreground_dat, LISA_combo):
    """Computes the LISA signal-to-noise ratio with inputs in SI units
    with a LISA mission of 4 yr and 2.5 million km arms
    Note that this does not take into account the Galactic binary foreground
    and assumes no orbital frequency evolution due to GW emission

    Parameters
    ----------
    m1 : Series
        primary mass [Msun]
    m2 : Series 
        secondary mass [Msun]
    porb : Series
        orbital period [s]
    ecc : Series
        eccentricity
    xGx, yGx, zGx : Series
        Galactocentric distance [kpc]
    n_harmonic : int
        number of frequency harmonics to include
    Tobs : float
        LISA observation time in seconds
    foreground_dat : DataFrame
        gw_freq [Hz], PSD [Hz^-1]   

    Returns
    -------
    SNR_dat : DataFrame
        DataFrame with columns ['freq', 'SNR'] where
        freq = gravitational-wave frequency in Hz and 
        SNR = LISA signal-to-noise ratio 
    """

    #LISA_combo_interp = LISA_combo(foreground_dat.freq) 
    dist = ((xGx-x_sun)**2+(yGx-y_sun)**2+(zGx-z_sun)**2)**0.5*parsec*1000.0
    mChirp = m_chirp(m1, m2)*Msun
    h_0 = 8*G/c**2*(mChirp)/(dist)*(G/c**3*2*np.pi*(1/(porb))*mChirp)**(2./3.)
    h_0_squared = h_0**2 
    if ecc.all() < 1e-4:
        n_harmonic = 2
        peters_g_factor = 0.5**2
        GW_freq_array = 2/porb
        LISA_curve_eval = LISA_combo(GW_freq_array)
        h_0_squared_dat = h_0_squared*peters_g_factor
        SNR = (h_0_squared_dat * Tobs / LISA_curve_eval)
    else:
        peters_g_factor = peters_gfac(ecc, n_harmonic)
        GW_freq_array = np.array([np.arange(1,n_harmonic)/p for p in porb]).T
        GW_freq_flat = GW_freq_array.flatten()
        LISA_curve_eval_flat = LISA_combo(GW_freq_flat)
        LISA_curve_eval = LISA_curve_eval_flat.reshape((n_harmonic-1,ecc.shape[0]))

        h_0_squared_dat = np.array([h*np.ones(n_harmonic-1) for h in h_0_squared]).T*peters_g_factor
        SNR_squared = h_0_squared_dat * Tobs / LISA_curve_eval
        SNR = (SNR_squared.sum(axis=0))**0.5
    gw_freq = peak_gw_freq(m1, m2, ecc, porb)    
    SNR_dat = pd.DataFrame(np.array([gw_freq, SNR]).T, columns=['gw_freq', 'SNR'])
    
    return SNR_dat


def LISA_PSD(m1, m2, porb, ecc, xGx, yGx, zGx, n_harmonic, Tobs):
    """Computes the LISA power spectral density with inputs in SI units
    with a LISA mission of Tobs [sec] and 2.5 million km arms
    Note that this does not take into account the Galactic binary foreground
    and assumes no orbital frequency evolution due to GW emission

    Parameters
    ----------
    m1 : Series 
        primary mass [msun]
    m2 : Series    
        secondary mass [msun]
    porb : Series
        orbital period [s]
    ecc : Series
        eccentricity
    xGx, yGx, zGx : Series
        Galactocentric distance [kpc]
    n_harmonic : int
        number of frequency harmonics to include
    Tobs : float
            LISA observation time [s]

    Returns
    -------
    PSD_dat : DataFrame
        DataFrame with columns ['freq', 'PSD'] where
        freq = gravitational wave frequency [Hz]
        PSD = LISA Power Spectral Density  
    """

    dist = ((xGx-x_sun)**2+(yGx-y_sun)**2+(zGx-z_sun)**2)**0.5*parsec*1000.0
    mChirp = m_chirp(m1, m2)*Msun

    h_0 = 8*G/c**2*(mChirp)/(dist)*(G/c**3*2*np.pi*(1/(porb))*mChirp)**(2./3.)
    h_0_squared = h_0**2
    if ecc.all() < 1e-4:
        n_harmonic = 2
        peters_g_factor = 0.5**2
        GW_freq_flat = 2/porb
        h_0_squared_dat = h_0_squared*peters_g_factor
        psd = h_0_squared_dat * Tobs
    else:
        peters_g_factor = peters_gfac(ecc, n_harmonic)
        GW_freq_array = np.array([np.arange(1,n_harmonic)/p for p in porb]).T
        GW_freq_flat = GW_freq_array.flatten()

        h_0_squared_dat = np.array([h*np.ones(n_harmonic-1) for h in h_0_squared]).T*peters_g_factor
        PSD = h_0_squared_dat * Tobs 
        psd = PSD.flatten()

    freq = GW_freq_flat
    PSD_dat = pd.DataFrame(np.vstack([freq, psd]).T, columns=['freq', 'PSD'])
    
    return PSD_dat 

def compute_foreground(psd_dat, Tobs=4*sec_in_year):
    """Computes the gravitational-wave foreground by binning the PSDs 
    in psd_dat according to the LISA frequency resolution where 1
    frequency bin has a binwidth of 1/(Tobs [s])

    Parameters
    ----------
    psd_dat : DataFrame
        DataFrame with columns ['freq', 'PSD'] where 
        freq = gravitational-wave frequency [Hz]
        PSD = LISA power spectral density
    Tobs (float):
        LISA observation time [s]; Default=4 yr

    Returns
    -------
    foreground_dat : DataFrame
        DataFrame with columns ['freq', 'PSD'] where 
        freq = gravitational-wave frequency of LISA frequency bins [Hz]
        PSD = LISA power spectral density of foreground
    """

    binwidth = 1.0/Tobs
    freqBinsLISA = np.arange(1e-6,1e-2,binwidth)
    binIndices = np.digitize(psd_dat.freq, freqBinsLISA)
    psd_dat['digits'] = binIndices
    power_sum = psd_dat[['PSD', 'digits']].groupby('digits').sum()['PSD']
    power_tot = np.zeros(len(freqBinsLISA)+1)
    power_tot[power_sum.index[:len(freqBinsLISA)-1]] = power_sum
    foreground_dat = pd.DataFrame(np.vstack([freqBinsLISA, power_tot[1:]]).T,\
                                  columns=['freq', 'PSD'])
    foreground_dat = foreground_dat.rolling(10).median()
    foreground_dat = foreground_dat.iloc[10:]
    return foreground_dat
