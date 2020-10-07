from scipy.interpolate import interp1d
from scipy.integrate import quad
from scipy.special import hyp2f1
from scipy.optimize import brentq
import numpy as np
from numpy.random import uniform, normal


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
        rOrvirMax /= 10.0

    def y_zero(r):
        return r / virial_radius_analytic(gamma, r) - r_max

    yMin = y_zero(rOrvirMin)
    yMax = y_zero(rOrvirMax)

    rmax_vir = brentq(y_zero, yMin, yMax)

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

    # Use an interpolator, but flip x and y for the CDF
    interpolator = interp1d(mass_enclosed_grid, radii_grid, kind="cubic")

    X = uniform(size=N)

    positions = interpolator(X)

    return positions


def get_velocities(r, r_max_cluster, gamma):
    """
    Uses the spherical Jeans functions to sample the velocity dispersion for the
    cluster at different radii, then draws a random, isotropic velocity for each
    star

    returns (vr,vt) with same length as r
    """

    N = len(r)

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

    vr = vx
    vt = np.sqrt(vy ** 2 + vz ** 2)

    return vr, vt


def draw_vr_vt_r(N=100000, r_max=300, gamma=4):
    """
    Draw random velocities and positions from the Elson profile.

    N     = number of stars
    r_max = maximum number of virial radii for the farthest star
    gamma = steepness function of the profile

    Note that gamma=4 is the Plummer profile

    returns (vr,vt,r) in G=M_cluster=1 units
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

    return vr, vt, r
