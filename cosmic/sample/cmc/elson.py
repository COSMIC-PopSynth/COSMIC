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

"""`Elson`
"""

from scipy.interpolate import interp1d, CubicSpline
from scipy.integrate import quad
from scipy.special import hyp2f1
from scipy.optimize import brentq
import numpy as np
from numpy.random import uniform, normal
from scipy.stats import maxwell

__author__ = "Carl Rodriguez <carllouisrodriguez@gmail.com>"
__credits__ = "Carl Rodriguez <carllouisrodriguez@gmail.com>"


def M_enclosed(r, gamma, rho_0):
    """
    Compute the mass enclosed in an Elson profile at radius r with slope
    gamma, central concentration rho_0, and assumed scale factor a = 1

    In practice, this is only used to sample the positions of stars, so rho_0 is
    just picked to normalize the distribution (i.e. rho_0 s.t. M_enclosed(rmax) = 1)
    """

    prefactor = 4 * np.pi * rho_0 * r * r * r / 3.0

    hypergeometric = hyp2f1(1.5, (gamma + 1.0) / 2.0, 2.5, -(r ** 2))

    return prefactor * hypergeometric

def phi_r(r, gamma, rho_0):
    """
    Compute the gravitational potential of an Elson profile at radius r with slope
    gamma, central concentration rho_0, and assumed scale factor a = 1

    Needed to compute the escape speed (for resampling stars)
    """

    prefactor = - 4 * np.pi * rho_0 / (gamma - 1)

    hypergeometric = hyp2f1(0.5, (gamma - 1.0) / 2.0, 1.5, -(r ** 2))

    return prefactor * hypergeometric

def rho_r(r, gamma, rho_0):
    """
    Compute the density of the Elson profile at radius r
    Best to use the same normalized rho_0 from M_enclosed
    """
    

    return rho_0 * pow(1 + r * r, -(gamma + 1.0) / 2.0)

def virial_radius_analytic(gamma, r_max):
    """
    Virial radius is best calculated directly, since rmax may be pretty far from
    infinity.  Directly integrate 4*pi*r*rho*m_enclosed from 0 to rMax to get the
    binding energy, then just divide 0.5 by that.
    """

    rho_0 = 1.0 / M_enclosed(r_max, gamma, 1)

    integral, error = quad(lambda r: rho_r(r, gamma, rho_0) * r * 4 * np.pi * M_enclosed(r, gamma, rho_0), 0, r_max)

    return 1 / 2.0 / integral


def find_rmax_vir(r_max, gamma):
    """
    This is a little tricky: because the virial radius of the Elson profile
    depends on the maximum radius, if the profile is very flat (e.g. gamma~2)
    then you need a large maximum radius to get a large number of virial radius
    in there.  This basically finds r such that  r / rVir(r) = r_max.  In other
    word, how far out do we need to go to have r_max number of virial radii in the
    sample.
    """
    rOrvirMin = 1.0
    rOrvirMax = 20000.0

    # rvir ~ 1 for gamma > 4, and the M_enclosed integral doesn't converge.
    # Here's a cheap fix for that...
    if gamma > 4:
        rOrvirMax /= 10


    def y_zero(r):
        return r / virial_radius_analytic(gamma, r) - r_max

    yMin = y_zero(rOrvirMin)
    yMax = y_zero(rOrvirMax)

    rmax_vir = brentq(y_zero, rOrvirMin,rOrvirMax)

    return rmax_vir


def find_sigma_sqr(r, r_max_cluster, gamma):
    """
    Find the 1D velocity dispersion at a given radius r using one of the
    spherial Jeans equations (and assuming velocity is isotropic)
    """

    rho_0 = 1.0 / M_enclosed(r_max_cluster, gamma, 1)

    jeans_integrand = (
        lambda rr: rho_r(rr, gamma, rho_0) * M_enclosed(rr, gamma, rho_0) / rr / rr
    )

    integral, error = quad(jeans_integrand, r, r_max_cluster)

    return integral / rho_r(r, gamma, rho_0)


def get_positions(N, r_max_cluster, gamma):
    """
    This one's easy: the mass enclosed function is just the CDF of the mass
    density, so just invert that, and you've got positions.
    """

    # First normalize the central density s.t. M_enc(r_max) = 1
    rho_0 = 1.0 / M_enclosed(r_max_cluster, gamma, 1)

    radii_grid = np.logspace(-3, np.log10(r_max_cluster), 1000)

    # Add r=0 to the array
    radii_grid[0] = 0

    mass_enclosed_grid = M_enclosed(radii_grid, gamma, rho_0)

    # if it's a steep profile, you'll get a lot of similar numbers (to machine precision)
    # remove them here
    mass_enclosed_grid, unique_idx = np.unique(mass_enclosed_grid,return_index=True)
    radii_grid = radii_grid[unique_idx]

    # Use an interpolator, but flip x and y for the CDF
    interpolator = interp1d(mass_enclosed_grid, radii_grid, kind="cubic")

    X = uniform(size=N)

    positions = interpolator(X)

    return positions

def get_velocities_old(r, r_max_cluster, gamma):
    """
    Uses the spherical Jeans functions to sample the velocity dispersion for the
    cluster at different radii, then draws a random, isotropic velocity for each
    star

    NOTE: this gives you correct velocity dispersion and produces a cluster in virial
    equilibrium, but it's not strictly speaking correct (the distribution should
    have some kurtosis, in addition to getting the variance of the Gaussian correct).
    You can see the disagreement in the tail of the distributions when sampling a
    Plummer sphere.  This is what mcluster does...

    Use the new get_velocities function, which generates a distribution function
    directly from rho and samples velocities from it

    This is kept around because it's super useful to sample from when the rejection sampling
    gets too slow at the last X samples (e.g. gamma ~ 10)

    returns (vr,vt) with same length as r
    """

    N = len(r)

    rho_0 = 1.0 / M_enclosed(r_max_cluster, gamma, 1)

    # Rather than calculate the integral for every stellar position, just
    # create an interpolator of 1000 or so points, then sample from the
    # interpolated curve
    radii_grid = np.logspace(np.log10(min(r) * 0.999), np.log10(max(r) * 1.001), 1000)
    sigma_sqr_grid = np.array(
        [find_sigma_sqr(rr, r_max_cluster, gamma) for rr in radii_grid]
    )

    interpolator = interp1d(radii_grid, sigma_sqr_grid, kind="cubic")

    # Draw the 1D velocity dispersions
    sigma = np.sqrt(interpolator(r))

    # Then sample an isotropic gaussian in vx,vy,vz using those dispersions
    vx = normal(scale=sigma, size=N)
    vy = normal(scale=sigma, size=N)
    vz = normal(scale=sigma, size=N)
    
    # Finally the rejection sampling: because of the Gaussian, some of the velocities
    # we draw will naturally be above the local escape speed of the cluster.  To avoid 
    # this, we flag them, and resample them from the velocity distribution
    while False:
        ve = np.sqrt(-2*phi_r(r, gamma, rho_0))
        escapers = (np.sqrt(vx**2 + vy**2 + vz**2) > ve)

        number_of_escapers = np.sum(escapers)
        if number_of_escapers == 0:
            break
        else:
            vx_temp = normal(scale=sigma[escapers], size=number_of_escapers)
            vy_temp = normal(scale=sigma[escapers], size=number_of_escapers)
            vz_temp = normal(scale=sigma[escapers], size=number_of_escapers)

        vx[escapers] = vx_temp
        vy[escapers] = vy_temp
        vz[escapers] = vz_temp

    return np.sqrt(vx**2 + vy**2 + vz**2) 

def get_velocities(r, r_max_cluster, gamma):
    """
    The correct way to generate velocities: generate the distribution function from rho,
    then use rejection sampling.

    returns (vr,vt) with same length as r
    """

    N = len(r)

    # If we actually want a plummer sphere, we can default to the simpler function
    # using the analytic distribution function
    if gamma == 4:
        return get_velocities_plummer(r,r_max_cluster)
    # Otherwise need to compute f(E) numerically

    # First get phi and rho over 100 points or so
    rho_0 = 1.0 / M_enclosed(max(r), gamma, 1)

    r_stensil = np.logspace(np.log10(min(r)),np.log10(max(r)),100)
    phi = phi_r(r_stensil,gamma,rho_0)
    rho = rho_r(r_stensil,gamma,rho_0)

    # Then construct a cubic spline of rho as a function of phi
    # Technically should be psi = -phi, but FITPACK only accepts increasing values
    rho_phi_spline = CubicSpline(phi,rho,extrapolate=False)

    f_E = []
    energies = []

    # Now integrate over energies to get the distribution function 
    for en in np.linspace(min(-phi),max(-phi),100):

        integrand = lambda psi: rho_phi_spline(-psi,2) / np.sqrt(en - psi)

        # Eqn 4.46b of Binney and Tremain (2nd Edition)
        # (4.140b in 1st edition)
        temp_f_E = ((quad(integrand,min(-phi)*1.01,en,limit=200)[0] +
                     rho_phi_spline(-min(-phi),1)/np.sqrt(en)) / 17.493418) # sqrt(8)*pi^2

        # Extrapolating derivatives from a cubic spline is hella dangerous
        # Better to have it throw NaNs and discard them
        if np.isnan(temp_f_E):
            continue
        else:
            f_E.append(temp_f_E)
            energies.append(en)

    # Set f(E) to zero at zero energy 
    f_E = np.array([0] + f_E)
    energies = np.array([0] + energies)

    # Make an interpolator for f(E) and compute f and phi at every star
    f_E_interp = interp1d(energies,f_E)
    phi_at_r = phi_r(r,gamma,rho_0)
    f_E_at_r = f_E_interp(-phi_at_r)

    # Make an initial guess for every velocity
    x_rand = uniform(size=N)
    y_rand = uniform(size=N)
    velocities = x_rand*np.sqrt(-2*phi_at_r)

    resample = np.ones(N,dtype=bool)
    number_to_resample = N

    # Now do the rejection sampling
    while number_to_resample > 0: 
        resample[resample] = (y_rand[resample]*f_E_at_r[resample] > 
                              x_rand[resample]**2 * f_E_interp(-phi_at_r[resample]*(1 - x_rand[resample]**2)))
        number_to_resample = np.sum(resample)
        x_rand[resample] = uniform(size=number_to_resample)
        y_rand[resample] = uniform(size=number_to_resample)

    velocities = x_rand * np.sqrt(-2*phi_at_r)

    theta = np.arccos(uniform(-1,1,N))

    vr = velocities*np.cos(theta)
    vt = velocities*np.sin(theta)

    return vr, vt

def get_velocities_plummer(r, r_max_cluster):
    """
    The correct way to generate velocities

    this is a special case for gamma=4, which is the plummer sphere (and has an analytic distribution function)

    returns (vr,vt) with same length as r
    """

    N = len(r)


    # Make an initial guess for every velocity
    x_rand = uniform(size=N)
    y_rand = uniform(size=N)
    vesc = np.sqrt(2/np.sqrt(1+r**2))
    velocities = x_rand * vesc

    resample = np.ones(N,dtype=bool)
    number_to_resample = N

    # Now do the rejection sampling; this is taken from Aarseth's original algorithm using the dist function
    while number_to_resample > 0: 
        resample[resample] = 0.1*y_rand[resample] > x_rand[resample]**2 * (1 - x_rand[resample]**2)**3.5
        number_to_resample = np.sum(resample)
        x_rand[resample] = uniform(size=number_to_resample)
        y_rand[resample] = uniform(size=number_to_resample)

    # Scale up to the local escape speed
    velocities = x_rand * vesc 

    theta = np.arccos(uniform(-1,1,N))

    vr = velocities*np.cos(theta)
    vt = velocities*np.sin(theta)

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


def draw_r_vr_vt(N=100000, r_max=300, gamma=4):
    """
    Draw random velocities and positions from the Elson profile.

    N     = number of stars
    r_max = maximum number of virial radii for the farthest star
    gamma = steepness function of the profile

    Note that gamma=4 is the Plummer profile

    returns (vr,vt,r) in G=M_cluster=1 units

    N.B. if you're getting a "Expect x to not have duplicates" error
    with a large gamma, use a smaller r_max. 
    """
    # First convert r_max into max number of  virial radii
    r_max = find_rmax_vir(r_max, gamma)

    # Then draw the positions from the cumulative mass function
    r = get_positions(N, r_max, gamma)

    # Sort the radii (needed later)
    r = np.sort(r)

    # Finally, draw the velocities using the radii and one of the Jeans
    # equations
    vr, vt = get_velocities(r, r_max, gamma)

    # And scale them to Henon units 
    r,vr,vt = scale_pos_and_vel(r,vr,vt)

    return r, vr, vt
