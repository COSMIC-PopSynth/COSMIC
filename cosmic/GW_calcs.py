# -*- coding: utf-8 -*-
# Copyright (C) Katelyn Breivik (2017 - 2019)
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
import scipy.integrate as integrate
from scipy.optimize import brentq
from lisa_sensitivity import lisa_characteristic_noise as LISA_hc
import pandas as pd
from scipy.interpolate import interp1d

__author__ = 'Katelyn Breivik <katie.breivik@gmail.com>'
__credits__ = 'Scott Coughlin <scott.coughlin@ligo.org>'
__all__ = ['sep_from_p', 'p_from_a', 'comoving_distance', 'luminosity_distance',
           'z_from_lum_distance', 'm_chirp', 'peak_gw_freq', 'peters_g', 'F_e',
           'n_max', 'hc2_circ', 'hc2', 'da_dt', 'de_dt', 'peters_evolution',
           'snr_calc', 'snr_chirping', 'LISA_foreground_combo', 'LISA_foreground']


G = 6.67384e-11
c = 2.99792458e8
parsec = 3.08567758e16
Rsun = 6.955e8
Msun = 1.9891e30
day = 86400.0
m_in_au = 1.496e+11
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

def sep_from_p(p,m1,m2):
    """Computes the separation from orbital period with KEPLER III

    Parameters
    ----------
    p : float/array
        orbital period [yr]
    m1 : float/array
        primary mass [kg]
    m2 : float/array
        secondary mass [kg]

    Returns
    -------
    sep : float/array
        separation [au]
    """

    sep_3 = p**2*(m1+m2)
    return sep_3**(1/3.)

def p_from_a(sep,m1,m2):
    """ Computes separation from orbital period with kepler III

    Parameters
    ----------
    sep : float/array
        separation [au]
    m1 : float/array
        primary mass [msun]
    m2 : float/array
        secondary mass [msun]

    Returns
    -------
    p : float/array
        orbital period [yr]
    """

    p_2 = sep**3/(m1+m2)
    return p_2**0.5

def comoving_distance(z):
    """ Computes the comoving distance given a redshift
    with cosmo params from `Bennet+2014 <http://adsabs.harvard.edu/abs/2014ApJ...794..135B>`_
    following `Ned Wright's site <http://www.astro.ucla.edu/~wright/CosmoCalc.html>`_

    Parameters
    ----------
    z : float
        redshift

    Returns
    -------
    D_co : float
       Comoving distance [Mpc]
    """

    h = 0.696
    omegaM = 0.286
    omegaK = 0.
    omegaL = 1 - 0.286
    dh = 3000. / h
    e = lambda zp: 1./np.sqrt(omegaM*(1+zp)**3 + omegaL)

    D_co = dh*integrate.quad(e,0,z)[0]
    return D_co

def luminosity_distance(z):
    """ Computes the luminosity distance given a redshift

    Parameters
    ----------
    z : float
        redshift

    Returns
    -------
    D_lum : float
        Luminosity distance [Mpc]
    """
    D_lum = (1+z)*comoving_distance(z)

    return D_lum

def z_from_lum_distance(d):
    """ Computes the redshift given a luminosity distance

    Parameters
    ----------
    d : float
        Luminosity distance [Mpc]

    Returns
    -------
    z_redshift : float
        Redshift
    """
    zero = lambda z: luminosity_distance(z) - d
    z_redshift = brentq(zero,0,5)
    return z_redshift


def m_chirp(m1, m2, **kwargs):
    """Computes the chirp mass in the units of mass supplied

    Parameters
    ----------
    m1 : float or array
        primary mass
    m2 : float or array
        secondary mass
    z : float or array
        redshift

    Returns
    -------
    mchirp : float or array
        chirp mass in units of mass supplied in the detector
        frame if redshift is provided, else redshift=0
    """
    z = 0
    for key, value in kwargs.items():
        if key == 'z': z = value
    m_chirp = (m1*m2)**(3./5.)/(m1+m2)**(1./5.)
    m_chirp_obs = m_chirp*(1+z)
    return m_chirp_obs

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

def peters_g(e, n):
    """Computes the factor g_n from Peters & Mathews 1963

    Parameters
    ----------
    e : float or array
        eccentricity
    n : int
        nth frequency harmonic

    Returns
    -------
    g_fac : array
        array of g_n
    """
    g_fac = (n**4 / 32.0)*((ss.jv((n-2), (n*e)) -\
                            2*e*ss.jv((n-1), (n*e)) +\
                            2.0/n*ss.jv((n), (n*e)) +\
                            2*e*ss.jv((n+1), (n*e)) -\
                            ss.jv((n+2), (n*e)))**2 +\
                           (1-e**2)*(ss.jv((n-2), (n*e)) -\
                                     2*ss.jv((n), (n*e)) +\
                                     ss.jv((n+2), (n*e)))**2 +\
                           (4/(3.0*n**2))*(ss.jv((n), (n*e)))**2)

    return g_fac

def F_e(e):
    """ Computes the Peters & Mathews 1963 F(e)
    factor

    Parameters
    ----------
    e : float or array
        eccentricity

    Returns
    -------
    F_e : array
        eccentricity factor
    """

    nominator = 1 + 73./24.*e**2 + 37./96.*e**4
    denominator = (1 - e**2)**(7./2.)
    F_e = nominator/denominator

    return F_e

def n_max(ecc):
    """ Computes the maximum number of harmonics
    needed to compute SNRs given an array of
    eccentricities.

    Parameters
    ----------
    ecc : array
        Eccentricity

    Returns
    -------
    n_max : array
        Maximum number of harmonics to use when
        computing eccentric LISA SNR
    """

    n_max = 2*np.ones(len(ecc))

    ind_05, = np.where((ecc > 0.001) & (ecc <=0.5))
    n_max[ind_05] = 20

    ind_07, = np.where((ecc > 0.5) & (ecc <=0.7))
    n_max[ind_07] = 50

    ind_08, = np.where((ecc > 0.7) & (ecc <=0.8))
    n_max[ind_08] = 150

    ind_085, = np.where((ecc > 0.8) & (ecc <= 0.85))
    n_max[ind_085] = 300

    ind_085_great, = np.where(ecc > 0.85)
    n_max[ind_085_great] = 500

    return n_max

def hc2_circ(m1, m2, f_orb, D):
    """ Computes the characterstic power (square of the
    characterstic strain) for a population of binaries
    in the Galaxy. This assumes the binaries are
    circular and redshift can be passed.

    Parameters
    ----------
    m1 : array/float
        primary mass [msun]
    m2 : array/float
        secondary mass [msun]
    f_orb : array/float
        orbital frequency [Hz]
    D : array/float
        Luminosity distance [Mpc]

    Returns
    -------
    hc2_circ : array/float
        characteristic power
    """
    # Set the redshift based on the luminosity distance
    z = z_from_lum_distance(D)

    # assumes SI units!
    f_n = 2*f_orb
    f_n_z = f_n/(1+z)
    m_c_z = m_chirp(m1, m2, z=z)*Msun
    front_fac = 2./3.*np.pi**(-4./3.)*G**(5./3.)/c**3
    variables = (m_c_z)**(5./3.)*(D*parsec*1e6)**(-2)*\
                (f_n_z)**(-1./3.)*(1+z)**(-2)
    hc2_circ = front_fac*variables
    return hc2_circ

def hc2(m1, m2, f_orb, D, e, n):
    """ Computes the characterstic power (square of the
    characterstic strain) for a population of binaries
    in the Galaxy. Redshift can be passed as a kwarg.

    Parameters
    ----------
    m1 : float
        primary mass [msun]
    m2 : float
        secondary mass [msun]
    f_orb : float
        orbital frequency [Hz]
    D : float
        Luminosity Distance [Mpc]
    e : float
        eccentricity
    n : int
        integer harmonic of orbital frequency
    Returns
    -------
    hc2_n : float
        characteristic power at the nth harmonic
    """

    # Set the redshift based on the luminosity distance
    z = z_from_lum_distance(D)

    # assumes SI units!
    f_n = n*f_orb
    f_n_z = f_n/(1+z)
    m_c_z = m_chirp(m1, m2, z=z)*Msun
    D = D * parsec * 1e6
    front_fac = 2./3.*np.pi**(-4./3.)*G**(5./3.)/c**3
    variables = (m_c_z)**(5./3.)*(D)**(-2)*(2./n)**(2./3)*\
                (f_n_z)**(-1./3.)*(1+z)**(-2)*peters_g(e,n)/F_e(e)
    hc2_n = front_fac*variables
    return hc2_n

def da_dt(m1, m2, a, e):
    """Computes da/dt from Peters and Matthews 1963;
    assumes all units are SI

    Parameters
    ----------
    m1 : float
        primary mass [kg]
    m2 : float
        secondary mass [kg]
    a : float
        semimajor axis [m]
    e : float
        eccentricity

    Returns
    -------
    da/dt : float
        semimajor axis rate of change [m/s]
    """

    front_fac = -64./5.*G**3/c**5
    nominator = m1*m2*(m1+m2)*F_e(e)
    denominator = a**3

    da_dt = front_fac*nominator/denominator
    return da_dt

def de_dt(m1, m2, a, e):
    """Computes de/dt from Peters and Matthews 1963;
    assumes all units are SI

    Parameters
    ----------
    m1 : float
        primary mass [kg]
    m2 : float
        secondary mass [kg]
    a : float
        semimajor axis [m]
    e : float
        eccentricity

    Returns
    -------
    de_dt : float
         eccentricity rate of change [1/s]
    """


    front_fac = -304./15.*G**3/c**5
    nominator = m1*m2*(m1+m2)*F_e(e)
    denominator = a**4

    de_dt = front_fac*nominator/denominator
    return de_dt


def peters_evolution(a_0,e_0,m1,m2,t_evol,nstep):
    """
    Computes the evolution of a binary with spacing according
    to the amount of time specified. Calculation is in SI;
    Returns values in units specified.

    Parameters
    ----------
    a_0 : float
        initial semimajor axis [au]
    e_0 : float
        initial eccentricity
    m1 : float
        primary mass [msun]
    m2 : float
        secondary mass [msun]
    t_evol : float
        evolution time [sec]
    nstep : int
        number of timesteps for evolution

    Returns
    -------
    a_array : array
        evolution of semimajor axis [au]
    e_array : array
        evolution of eccentricity
    t_array : array
        evolution of time [yr]
    """
    t_start = 0.0
    t_end = t_evol
    a_0 = a_0*m_in_au
    m1 = m1*Msun
    m2 = m2*Msun
    a_array = []
    e_array = []
    t_array = []

    times = np.linspace(t_start, t_end, int(nstep))
    deltas = times[1:] - times[:-1]

    for t, delta in zip(times[1:], deltas):
        a_array.append(a_0/m_in_au)
        e_array.append(e_0)
        t_array.append(t/sec_in_year)
        if (e_0 > 0):
            a = a_0 + da_dt(m1, m2, a_0, e_0)*delta
            e = e_0 + de_dt(m1, m2, a_0, e_0)*delta
        elif a_0 > 0:
            a = a_0 + da_dt_circ(m1, m2, a_0)*delta
            e = e_0
        else:
            break

        a_0 = a
        e_0 = e
    return np.array(a_array), np.array(e_array), np.array(t_array)

def snr_calc(f_vals, hc_n_vals, noise_vals, averaging_factor=4./5.):
    """Computes the SNR sum over the characteristic strain evolution

    Parameters
    ----------
    f_vals : array
        frequencies in the characteristic strain evolution
    hc_n_vals : array
        characterstic strain evolution at different Peters harmonics
    noise_vals : array
        LISA sensitivity curve values
    averaging_factor : float
        Factor to average over inclination values of the binary

    Returns
    -------
    snr : float
        Sum of SNRs across all Peters harmonics
    """

    snr_squared_per_harm = averaging_factor*np.trapz(1./f_vals*(hc_n_vals/noise_vals)**2., x=f_vals)
    snr = np.sqrt(snr_squared_per_harm.sum())
    return snr


def snr_chirping(m1, m2, a, e, d_lum, t_obs):
    """Computes the LISA SNR for a chirping, eccentric binary
    using the Robson, Cornish, and Liu 2018 characteristic
    strain sensivity curve

    Parameters
    ----------
    m1 : float
        primary mass [msun]
    m2 : float
        secondary mass [msun]
    a : float
        semimajor axis [au]
    e : float
        eccentricity
    d_lum : float
        Distance [Mpc]
    t_obs : float
        LISA observation time [yr]

    Returns
    -------
    snr : float
        LISA snr
    """
    LISA_interp = LISA_hc()
    z = z_from_lum_distance(d_lum)
    t_obs = t_obs * sec_in_year
    a_evol, e_evol, t_evol = peters_evolution(a,e,m1,m2,t_obs,1e3)
    porb = p_from_a(a, m1, m2) * sec_in_year
    f_orb_start = 1/porb
    f_orb_evol = 1/(p_from_a(a_evol, m1, m2) * sec_in_year)
    t_evol_log = np.logspace(-6, np.log10(t_obs), 5000)
    forb_evol_log = np.interp(t_evol_log,
                              xp = t_evol * sec_in_year,
                              fp = f_orb_evol)
    e_evol_log = np.interp(t_evol_log,
                           xp = t_evol * sec_in_year,
                           fp = e_evol)
    n_harm = int(n_max(np.array([e]))[0])
    h_c_squared = []
    freqs = []
    for n in range(1,n_harm):
        h_c_squared.append(hc2(m1, m2, forb_evol_log, d_lum, e_evol_log, n))
        freqs.append(forb_evol_log*n)

    h_c_squared = np.array(h_c_squared)
    freqs = np.array(freqs)
    snr = snr_calc(freqs*(1+z), h_c_squared**0.5, LISA_interp(freqs*(1+z)), 1)

    return snr

def LISA_foreground_combo(foreground_freq, foreground_hc2):
    """Computes the combination of the population foreground and LISA
    sensitivity in characteristic strain

    Parameters
    ----------
    foreground_freq : array
        gravitational wave frequencies [Hz]
    foreground_hc2 : array
        characteristic powers [Hz^-1]

    Returns
    -------
    LISA_curve_combo : interpolation
        interpolation of characteristic power for the Galactic
        foreground and LISA sensitivity combined
    """
    LISA_interp = LISA_hc()

    LISA_at_foreground = LISA_interp(foreground_freq)**2
    power_combo = foreground_hc2 + LISA_at_foreground

    above_foreground_freq = np.logspace(np.log10(max(foreground_freq)+1e-4, 1e-1, 100))
    LISA_above_foreground = LISA_interp(above_foreground_freq)**2

    freqs_total = np.hstack([foreground_freq, above_foreground_freq])
    hc2_total = np.hstact([power_combo, LISA_above_foreground])

    LISA_combo_interp = interp1d(freqs_total, hc_total)

    return LISA_combo_interp

def LISA_foreground(m1, m2, f_orb, d_lum, t_obs, rolling_window=100):
    """Computes the gravitational-wave foreground in characteristic
    power according to the LISA frequency resolution where 1
    frequency bin has a binwidth of 1/(Tobs [s])

    Parameters
    ----------
    m1 : array
        primary mass [msun]
    m2 : array
        secondary mass [msun]
    f_orb : array/float
        orbital frequency [Hz]
    d_lum : array
        Distance [Mpc]
    z : array
        redshift, Tobs=4):
    t_obs : float
        LISA observation time [yr]; Default=4 yr
    rolling_window : int
        Number of bins to compute a rolling median over

    Returns
    -------
    foreground_freq : array
        array of frequencies of the LISA foreground
    foreground_hc2 : array
        array of characteristic power of the LISA foreground
    """
    t_obs = t_obs * sec_in_year

    binwidth = 1.0/t_obs
    freq_bins_LISA = np.arange(1e-8,5e-3,binwidth)
    z_redshift = zAtLuminosityDistance(d_lum)
    hc2_array = hc2_circ(m1, m2, f_orb, d_lum, z=z_redshift)
    freq_array = 2*f_orb
    bin_indices = np.digitize(freq_array, freq_bins_LISA)
    psd_dat = pd.DataFrame(np.vstack([freq_array, hc_array, bin_indices]),\
                           columns=['f_gw', 'hc2', 'bin_digits'])

    power_sum = psd_dat[['hc2', 'bin_digits']].groupby('bin_digits').sum()['hc2']
    power_tot = np.zeros(len(freq_bins_LISA)+1)
    power_tot[power_sum.index[:len(freq_bins_LISA)-1]] = power_sum
    foreground_dat = pd.DataFrame(np.vstack([freq_bins_LISA, power_tot[1:]]).T,\
                                  columns=['freq', 'hc2'])
    foreground_dat = foreground_dat.rolling(rolling_window).median()
    foreground_dat = foreground_dat.iloc[rolling_window:]

    return np.array(foreground_dat['freq']), np.array(foreground_dat['hc2'])
