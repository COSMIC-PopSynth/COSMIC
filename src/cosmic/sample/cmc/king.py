# -*- coding: utf-8 -*-
# Copyright (C) Carl Rodriguez (2020 - 2021)
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

"""`King`
"""

import numpy as np
from numpy.random import uniform, normal
from scipy.integrate import RK45, quad, simpson, cumulative_trapezoid
from scipy.interpolate import interp1d
from scipy.special import erf

__author__ = "Carl Rodriguez <carllouisrodriguez@gmail.com>"
__credits__ = "Carl Rodriguez <carllouisrodriguez@gmail.com>"
__all__ = ["calc_rho", "integrate_king_profile", "virial_radius_numerical", "find_sigma_sqr", "get_positions", "get_velocities", "scale_pos_and_vel", "draw_r_vr_vt"]

def calc_rho(w):
    """
    Returns the density (unnormalized) given w = psi/sigma^2
    Note that w(0) = w_0, the main parameter of a King profile
    """
    rho = np.exp(w)*erf(np.sqrt(w)) - np.sqrt(4.0*w/np.pi)*(1.0+2.0*w/3.0)
    return rho

def integrate_king_profile(w0,tidal_boundary=1e-6):
    """ 
    Integrate a King Profile of a given w_0 until the density (or phi/sigma^2) 
    drops below tidal_boundary limit (1e-8 times the central density by default)
    
    Let's define some things:
    The King potential is often expressed in terms of w = psi / sigma^2 (psi is phi0 - phi, so just positive potential)
    (note that sigma is the central velocity dispersion with an infinitely deep potential, and close otherwise)
    
    in the center, w_0 = psi_0 / sigma^2 (the free parameter of the King profile)
    
    The core radius is defined (King 1966) as r_c = sqrt(9 * sigma^2 / (4 pi G rho_0))
    
    If we define new scaled quantities
            r_tilda = r/r_c
            rho_tilda = rho/rho_o
    We can rewrite Poisson's equation, (1/r^2) d/dr (r^2 dphi/dr)  =  4 pi G rho as:
            d^2(r_tilda w)/dr_tilda^2 = 9 r_tilda rho_tilda
    
    After that, all we need is initial conditions:
    w(0) = w_0
    w'(0) = 0

    returns (radii, rho, phi, M_enclosed)
    """

    rho_0 = calc_rho(w0)

    def ode_rhs(r_tilda,w_vec):
        ## Unpack the w_vector (note I'm cheating by putting rho_0 as an element)
        w,w_dot,rho_0 = w_vec

        ## Compute rho and second derivative of w
        rho = calc_rho(w)
        w_ddot = -9*rho/rho_0 - 2*w_dot/r_tilda

        ## return dw/dr_tilda, d^2w/dr_tilda^2, d(rho_tilda)/dr_tilda = 0
        return (w_dot,w_ddot,0)

    king_profile = RK45(ode_rhs,1e-4,[w0,0,rho_0],1e5,rtol=1e-10,atol=1e-12)

    at_the_tidal_boundary = False
    r = []
    rho_r = []
    phi_r = []
    
    while not at_the_tidal_boundary:
        w,w_dot,rho_0 = king_profile.y
        rho = calc_rho(w)

        r.append(king_profile.t)
        rho_r.append(rho)
        phi_r.append(w)

        if rho/rho_0 < tidal_boundary:
            at_the_tidal_boundary = True
        else:
            king_profile.step()
            
    r = np.array(r)
    rho_r = np.array(rho_r)
    
    ## Finally compute the cumulative mass enclosed
    M_enclosed = cumulative_trapezoid(4*np.pi*r**2*rho_r, r, initial=0)
    
    return r, rho_r, phi_r, M_enclosed


def virial_radius_numerical(r, rho_r, M_enclosed):
    """
    Virial radius is best calculated directly.  Directly integrate 
    4*pi*r*rho*m_enclosed over the samples binding energy, then just divide 0.5 by that.
    """
    integral = simpson(y=rho_r * r * 4 * np.pi * M_enclosed,x=r)

    return 1 / 2.0 / integral

def find_sigma_sqr(r_sample, r, rho_r, M_enclosed):
    """
    Find the 1D velocity dispersion at a given radius r using one of the
    spherial Jeans equations (and assuming velocity is isotropic)
    """

    ## This needs to be done continuously, so use interpolators
    rho_interp = interp1d(r,rho_r)
    M_enc_interp = interp1d(r,M_enclosed)
    
    jeans_integrand = lambda rr: rho_interp(rr) * M_enc_interp(rr) / rr / rr
    jeans_integrand = interp1d(r,rho_r*M_enclosed/r/r)

    integral, error = quad(jeans_integrand, r_sample, r[-1])
    return integral / rho_interp(r_sample)

def get_positions(N, r, M_enclosed):
    """
    This one's easy: the mass enclosed function is just the CDF of the mass
    density, so just invert that, and you've got positions.

    Note that here we've already normalized M_enclosed to 1 at r_tidal
    """

    # Use an interpolator, but flip x and y for the CDF
    interpolator = interp1d(M_enclosed, r, kind="cubic")

    X = uniform(size=N)

    positions = interpolator(X)

    return positions

def get_velocities(r, r_profile, psi_profile, M_enclosed_profile):
    """
    The correct way to generate velocities: start from the distribution function
    and use rejection sampling.

    returns (vr,vt) with same length as r
    """

    N = len(r)

    # create an interpolator for the psi
    psi_r = interp1d(r_profile,psi_profile)
    psi = psi_r(r)

    # Make an initial guess for every velocity (some fraction of the escape speed)
    x_rand = uniform(size=N)
    y_rand = uniform(size=N)

    resample = np.ones(N,dtype=bool)
    number_to_resample = N

    minF = (np.exp(psi) - 1)*2*psi
    v = x_rand * np.sqrt(2*psi)
    f_0 = y_rand*minF
    f = (np.exp(psi-v**2/2)-1)*v**2

    # Now do the rejection sampling
    while number_to_resample > 0: 
        v[resample] = x_rand[resample] * np.sqrt(2*psi[resample])
        f_0[resample] = y_rand[resample]*minF[resample]
        f[resample] = (np.exp(psi[resample]-v[resample]**2/2)-1)*v[resample]**2

        resample[resample] = f[resample] < f_0[resample]
        number_to_resample = np.sum(resample)
        x_rand[resample] = uniform(size=number_to_resample)
        y_rand[resample] = uniform(size=number_to_resample)

    theta = np.arccos(uniform(-1,1,N))

    vr = v*np.cos(theta)
    vt = v*np.sin(theta)

    return vr, vt


def scale_pos_and_vel(r,vr,vt):
    """
    Scale the positions and velocities to be in N-body units
    If we add binaries we'll do this again in initialcmctable.py

    takes r, vr, and vt as input

    returns (r,vr,vt) scaled to Henon units
    """

    # Normalize masses to 1/Mtotal
    mass = np.ones_like(r)/len(r)
    cumul_mass = np.cumsum(mass)

    radius = np.array(r)
    radius_p1 = np.append(radius[1:], [1e100])

    # Then compute the total kinetic and potential energy
    # There's probably a cleaner way to do the PE (this is a one-line version
    #  of the for loop we use in CMC; vectorized and pythonic, but sloppy)
    KE = 0.5 * sum(mass * (vr ** 2 + vt ** 2))
    PE = 0.5 * sum(
        mass[::-1]
        * np.cumsum((cumul_mass * (1.0 / radius - 1.0 / radius_p1))[::-1])
    )

    # Compute the position and velocity scalings
    rfac = 2 * PE
    vfac = 1.0 / np.sqrt(4 * KE)

    return r*rfac, vr*vfac, vt*vfac

def draw_r_vr_vt(N=100000, w_0=6, tidal_boundary=1e-6):
    """
    Draw random velocities and positions from the King profile.

    N              = number of stars
    w_0            = King concentration parameter (-psi/sigma^2) 
    tidal_boundary = ratio of rho/rho_0 where we truncate the tidal boundary

    returns (vr,vt,r) in G=M_cluster=1 units
    """
    # First integrate the differential equation of a King profile 
    r_profile, rho_profile, psi_profile, M_enclosed_profile = integrate_king_profile(w_0, tidal_boundary=tidal_boundary)

    ## Normalize the masses to 1 at r_tidal
    rho_profile /= M_enclosed_profile[-1]
    M_enclosed_profile /= M_enclosed_profile[-1]

    ## Normalize r (currently in units of the King core radius) to the virial radius
    r_profile /= virial_radius_numerical(r_profile, rho_profile, M_enclosed_profile)

    # Then draw the positions from the cumulative mass function
    r = get_positions(N, r_profile, M_enclosed_profile)

    # Sort the radii (needed later)
    r = np.sort(r)

    # Finally, draw the velocities using the radii and one of the Jeans
    # equations
    vr, vt = get_velocities(r, r_profile, psi_profile, M_enclosed_profile)

    # Scale into Henon units
    r, vr, vt = scale_pos_and_vel(r,vr,vt)

    return r, vr, vt
