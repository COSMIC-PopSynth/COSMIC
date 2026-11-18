"""Unit test for cosmic kick routine
"""

__author__ = 'Tom Wagg <tomjwagg@gmail.com>'

import os
import unittest
import numpy as np
from scipy.stats import maxwell, norm
from scipy.optimize import brentq
import pandas as pd

from cosmic.sample.initialbinarytable import InitialBinaryTable, INITIAL_CONDITIONS_COLUMNS_ALL
from cosmic.evolve import Evolve, INITIAL_BINARY_TABLE_SAVE_COLUMNS, FLATTENED_NATAL_KICK_COLUMNS

import warnings
warnings.filterwarnings("ignore")


TEST_DATA_DIR = os.path.join(os.path.split(__file__)[0], 'data')
INIT_CONDITIONS = pd.read_hdf(os.path.join(TEST_DATA_DIR, 'initial_conditions_for_testing.hdf5'), key='initC')
SSEDict = {"stellar_engine": "sse", "path_to_tracks": "", "path_to_he_tracks": ""}
BSEFlag_columns = list(set(INITIAL_BINARY_TABLE_SAVE_COLUMNS) - set(INITIAL_CONDITIONS_COLUMNS_ALL))
BSEDict = INIT_CONDITIONS[BSEFlag_columns].to_dict(orient='index')[0]
BSEDict['qcrit_array'] = [0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]
BSEDict['natal_kick_array'] = [[-100.0,-100.0,-100.0,-100.0,0.0], [-100.0,-100.0,-100.0,-100.0,0.0]]
BSEDict['fprimc_array'] = [2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,
                           2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,
                           2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0]

# avoid using the same randomseed
del BSEDict["bin_num"], BSEDict["randomseed"]


def evolve_with_first_sn_kicks(binary, natal_kick, phi, theta, mean_anomaly, tphysf):
    """Evolve copies of a binary, each with its own supplied first-SN kick

    Each argument (other than ``binary`` and ``tphysf``) is an array with one
    entry per copy, in the units of ``natal_kick_array``. All copies use the
    same random seed, so they reach the SN in an identical state.
    """
    n = len(mean_anomaly)
    binaries = binary.loc[binary.index.repeat(n)].reset_index(drop=True)
    binaries["tphysf"] = tphysf
    for col, value in [("natal_kick", natal_kick), ("phi", phi), ("theta", theta),
                       ("mean_anomaly", mean_anomaly), ("randomseed", 0.0)]:
        binaries[col + "_1"] = value
        binaries[col + "_2"] = 0.0 if col == "randomseed" else -100.0
    # without natal_kick_array (or its flattened columns) in the BSEDict, the
    # per-binary columns are used
    bse = {k: v for k, v in BSEDict.items()
           if k != "natal_kick_array" and k not in FLATTENED_NATAL_KICK_COLUMNS}
    bse["kickflag"] = 1
    bpp, _, _, kick_info = Evolve.evolve(initialbinarytable=binaries, BSEDict=bse,
                                         SSEDict=SSEDict, nproc=1, randomseed=42)
    return bpp, kick_info[kick_info["star"] == 1]


def wrap_degrees(angle):
    """Wrap an angle (in degrees) to [-180, 180)"""
    return (np.asarray(angle) + 180.0) % 360.0 - 180.0


def kepler_state(a, e, mean_anomaly, gm):
    """Relative position and velocity in a Keplerian orbit

    Periapsis is along +x and the motion is along +y at periapsis, matching
    the frame that kick.f uses. ``mean_anomaly`` is in radians.
    """
    E = brentq(lambda E: E - e * np.sin(E) - mean_anomaly, 0.0, 2 * np.pi,
               xtol=1e-15, rtol=1e-15)
    n = np.sqrt(gm / a**3)
    r = a * (1 - e * np.cos(E))
    pos = np.array([a * (np.cos(E) - e), a * np.sqrt(1 - e**2) * np.sin(E), 0.0])
    vel = n * a**2 / r * np.array([-np.sin(E), np.sqrt(1 - e**2) * np.cos(E), 0.0])
    return pos, vel


def euler_rotation(theta, phi, psi):
    """z-x-z rotation matrix R = Rz(phi) Rx(theta) Rz(psi), angles in degrees

    This matches ChangeBasis in kick.f, which maps a vector from the new
    frame to the old one.
    """
    theta, phi, psi = np.radians([theta, phi, psi])

    def rot_z(a):
        return np.array([[np.cos(a), -np.sin(a), 0], [np.sin(a), np.cos(a), 0], [0, 0, 1]])

    rot_x = np.array([[1, 0, 0], [0, np.cos(theta), -np.sin(theta)],
                      [0, np.sin(theta), np.cos(theta)]])
    return rot_z(phi) @ rot_x @ rot_z(psi)


class TestKick(unittest.TestCase):
    """`TestCase` for the cosmic kick routine
    """

    def test_disberg(self):
        """Test the Disberg kick routine
        """
        # define a simple single star that will result in a kick
        single_star = InitialBinaryTable.InitialBinaries(
            m1=20 + np.random.rand() * 0.01, m2=0.0, porb=0,
            ecc=-1, tphysf=100.0, kstar1=1, kstar2=15,
            metallicity=0.02
        )
        N = 10000
        ibt = single_star.loc[single_star.index.repeat(N)].reset_index()

        # turn off ECSN and bh fallback
        BSEDict["ecsn"] = 0
        BSEDict["ecsn_mlow"] = 0
        BSEDict["bhflag"] = 3

        # evolve using disberg
        BSEDict["kickflag"] = 5
        _, _, _, kick_info = Evolve.evolve(initialbinarytable=ibt, BSEDict=BSEDict, SSEDict=SSEDict, nproc=1)
        natal_kicks_disberg = kick_info['natal_kick'][kick_info['natal_kick'] != 0.0]

        # fit a lognormal distribution and ensure it matches the expected values
        # (mu = 5.60, sigma = 0.68) to within 5 standard errors of the fit:
        # SE(mu) = sigma / sqrt(n) and SE(sigma) ~ sigma / sqrt(2n)
        mu_d, sigma_d = norm.fit(np.log(natal_kicks_disberg))
        n = len(natal_kicks_disberg)
        self.assertGreater(n, 0.9 * N)
        self.assertAlmostEqual(mu_d, 5.60, delta=5 * 0.68 / np.sqrt(n))
        self.assertAlmostEqual(sigma_d, 0.68, delta=5 * 0.68 / np.sqrt(2 * n))

    def test_hobbs(self):
        """Test the Hobbs kick routine
        """
        # define a simple single star that will result in a kick
        single_star = InitialBinaryTable.InitialBinaries(
            m1=20 + np.random.rand() * 0.01, m2=0.0, porb=0,
            ecc=-1, tphysf=100.0, kstar1=1, kstar2=15,
            metallicity=0.02
        )
        N = 10000
        ibt = single_star.loc[single_star.index.repeat(N)].reset_index()

        # turn off ECSN and bh fallback
        BSEDict["ecsn"] = 0
        BSEDict["ecsn_mlow"] = 0
        BSEDict["bhflag"] = 3

        # evolve using hobbs
        BSEDict["kickflag"] = 1
        _, _, _, kick_info = Evolve.evolve(initialbinarytable=ibt, BSEDict=BSEDict, SSEDict=SSEDict, nproc=1)
        natal_kicks = kick_info['natal_kick'][kick_info['natal_kick'] != 0.0]

        # fit a maxwellian to the hobbs natal kicks and ensure it matches the
        # expected dispersion (sigma = 265 km/s) to within 5 standard errors of
        # the fit: SE(sigma) = sigma / sqrt(6n)
        _, s_hobbs = maxwell.fit(natal_kicks, floc=0.0)
        n = len(natal_kicks)
        self.assertGreater(n, 0.9 * N)
        self.assertAlmostEqual(s_hobbs, 265.0, delta=5 * 265.0 / np.sqrt(6 * n))

    def test_kick_collision(self):
        """Test that a natal kick aimed at the companion destroys it

        A He star + MS binary goes through a SN from a circular orbit. We read
        the pre-SN state from the bpp, then re-evolve with a natal kick chosen
        so that the post-SN relative velocity points almost straight at the
        companion: the stars are approaching and the pericentre of the
        (unbound) post-SN orbit is ~0.5 Rsun, well inside the companion's
        ~4 Rsun radius, so this must end in a collision. The mirrored kick,
        which sends the remnant directly away from the companion, must instead
        disrupt the binary and leave the companion intact.
        """
        G = 1.3271244e11  # km^3 / (Msun s^2), as in kick.f
        rsunkm = 6.96e5
        binary = InitialBinaryTable.InitialBinaries(
            m1=6.0, m2=10.0, porb=5.0, ecc=0.0, tphysf=2.0,
            kstar1=7, kstar2=1, metallicity=0.02
        )
        bse = dict(BSEDict, kickflag=1)

        # first pass: get the state of the binary right before the SN
        bpp, _, _, _ = Evolve.evolve(initialbinarytable=binary, BSEDict=bse,
                                     SSEDict=SSEDict, nproc=1, randomseed=42)
        pre_sn = bpp[bpp["evol_type"] == 15].iloc[0]
        self.assertEqual(pre_sn["kstar_2"], 1)
        self.assertGreater(pre_sn["sep"], 0.0)
        self.assertEqual(pre_sn["ecc"], 0.0)

        # with a mean anomaly of 0, the exploding star sits at (a, 0, 0)
        # relative to its companion and moves with velocity (0, v_orb, 0)
        a = pre_sn["sep"] * rsunkm
        v_orb = np.sqrt(G * (pre_sn["mass_1"] + pre_sn["mass_2"]) / a)

        # target post-SN relative velocities (in-plane, km/s)
        for v_rel, expect_collision in [((-1500.0, 50.0), True),
                                        ((1500.0, 50.0), False)]:
            with self.subTest(v_rel=v_rel):
                kick = np.array([v_rel[0], v_rel[1] - v_orb])
                natal_kick_array = [
                    [np.linalg.norm(kick), 0.0,
                     np.degrees(np.arctan2(kick[1], kick[0])) % 360.0, 0.0, 0.0],
                    [-100.0, -100.0, -100.0, -100.0, 0.0]
                ]
                bpp, _, _, kick_info = Evolve.evolve(
                    initialbinarytable=binary,
                    BSEDict=dict(bse, natal_kick_array=natal_kick_array),
                    SSEDict=SSEDict, nproc=1, randomseed=42
                )
                post_sn = bpp.iloc[bpp["evol_type"].tolist().index(15) + 1]
                self.assertEqual(kick_info["disrupted"].iloc[0], 1)
                if expect_collision:
                    # coalescence, companion destroyed
                    self.assertEqual(post_sn["evol_type"], 6)
                    self.assertEqual(post_sn["kstar_2"], 15)
                else:
                    # disrupted, companion survives
                    self.assertEqual(post_sn["evol_type"], 11)
                    self.assertEqual(post_sn["kstar_2"], 1)

    def test_kepler_high_eccentricity(self):
        """Test the eccentric anomaly solve for highly eccentric orbits

        Plain Newton-Raphson on Kepler's equation (starting from E = M) can
        cycle forever for eccentricities close to 1. Here a He star + MS
        binary goes SN with e ~ 0.9999 and *no* natal kick, at a range of
        supplied mean anomalies. These include a few where the old iteration
        never converged for this eccentricity (at least with gfortran on
        arm64 macOS; exactly which mean anomalies cycle depends on rounding),
        and some either side of the edge of the "critical region" near
        periapsis (M < 0.0045 rad ~ 0.258 deg from it), inside which the
        solver switches to bisection.

        With no natal kick, the change in systemic velocity is purely the
        mass-loss term, v_cm ∝ -v_rel ∝ (sin E, -sqrt(1 - e^2) cos E, 0), so
        we can recover the eccentric anomaly E from kick_info and check that
        it satisfies Kepler's equation for the supplied mean anomaly.
        """
        binary = InitialBinaryTable.InitialBinaries(
            m1=6.0, m2=10.0, porb=4.8e6, ecc=0.9999, tphysf=2.0,
            kstar1=7, kstar2=1, metallicity=0.02
        )
        critical = np.array([1e-4, 0.01, 0.1, 0.25, 0.26, 0.3, 1.0])
        mean_anomalies = np.concatenate([np.arange(0.0, 360.1, 5.0),
                                         [11.74, 18.87, 20.24, 332.91, 341.18],
                                         critical, 360.0 - critical])
        n = len(mean_anomalies)
        bpp, kick_info = evolve_with_first_sn_kicks(
            binary, natal_kick=np.zeros(n), phi=np.zeros(n), theta=np.zeros(n),
            mean_anomaly=mean_anomalies, tphysf=1.2
        )

        eccs, ecc_anoms = [], []
        for bin_num, ki in kick_info.iterrows():
            this_bpp = bpp.loc[[bin_num]].reset_index(drop=True)
            sn_row = this_bpp["evol_type"].tolist().index(15)
            ecc = this_bpp.loc[sn_row, "ecc"]
            self.assertGreater(ecc, 0.9998)
            eccs.append(ecc)

            v_1 = ki[["delta_vsysx_1", "delta_vsysy_1", "delta_vsysz_1"]].values
            if ki["disrupted"] == 1:
                # each star's velocity is given, combine to get the CoM velocity
                m_1 = this_bpp.loc[sn_row + 1, "mass_1"]
                m_2 = this_bpp.loc[sn_row, "mass_2"]
                v_2 = ki[["delta_vsysx_2", "delta_vsysy_2", "delta_vsysz_2"]].values
                v_cm = (m_1 * v_1 + m_2 * v_2) / (m_1 + m_2)
            else:
                v_cm = v_1
            ecc_anoms.append(np.arctan2(v_cm[0], -v_cm[1] / np.sqrt(1 - ecc**2)))
            self.assertAlmostEqual(v_cm[2], 0.0, places=8)

        eccs, ecc_anoms = np.array(eccs), np.array(ecc_anoms)
        mean_anom_check = np.degrees(ecc_anoms - eccs * np.sin(ecc_anoms))
        # (the solver itself is accurate to ~1e-13 deg, reconstructing E from
        # the velocities limits this check to ~1e-11 deg)
        np.testing.assert_allclose(wrap_degrees(mean_anom_check - kick_info["mean_anomaly"]),
                                   0.0, atol=1e-9)
        np.testing.assert_allclose(kick_info["mean_anomaly"], mean_anomalies)

    def test_euler_angles_circular(self):
        """Test the orbital-plane Euler angles after a SN in a circular orbit

        In a circular pre-SN orbit, the separation at mean anomaly M is
        r = a (cos M, sin M, 0) in the frame used for the kick_info
        velocities. A kick with an out-of-plane component tilts the orbit
        about r, so the line of nodes (z x h) points along +r for an upward
        kick and -r for a downward one: phi_euler = M or M - 180 deg. The tilt
        is the angle between h = r x v_rel and the z-axis.

        These angles are needed to combine velocities from the two SNe, and
        previously phi_euler was drawn at random for any circular orbit.
        """
        G = 1.3271244e11  # km^3 / (Msun s^2), as in kick.f
        rsunkm = 6.96e5
        binary = InitialBinaryTable.InitialBinaries(
            m1=6.0, m2=10.0, porb=5.0, ecc=0.0, tphysf=2.0,
            kstar1=7, kstar2=1, metallicity=0.02
        )
        mean_anomalies = np.array([0.0, 37.0, 60.0, 135.0, 200.0, 290.0, 330.0])
        kick_phis = np.array([30.0, -30.0, 60.0, -10.0, 45.0, 80.0, -75.0])
        kick_thetas = np.array([0.0, 90.0, 200.0, 45.0, 300.0, 120.0, 10.0])
        n = len(mean_anomalies)
        bpp, kick_info = evolve_with_first_sn_kicks(
            binary, natal_kick=np.full(n, 50.0), phi=kick_phis, theta=kick_thetas,
            mean_anomaly=mean_anomalies, tphysf=1.2
        )
        # all should survive the SN, otherwise the Euler angles are random
        self.assertTrue((kick_info["disrupted"] == 0).all())

        pre_sn = bpp[bpp["evol_type"] == 15]
        self.assertTrue((pre_sn["ecc"] == 0.0).all())
        a = pre_sn["sep"].values * rsunkm
        v_orb = np.sqrt(G * (pre_sn["mass_1"].values + pre_sn["mass_2"].values) / a)

        M = np.radians(mean_anomalies)
        k_phi, k_theta = np.radians(kick_phis), np.radians(kick_thetas)
        r = a[:, None] * np.stack([np.cos(M), np.sin(M), np.zeros(n)], axis=1)
        v_rel = (v_orb[:, None] * np.stack([-np.sin(M), np.cos(M), np.zeros(n)], axis=1)
                 + 50.0 * np.stack([np.cos(k_phi) * np.cos(k_theta),
                                     np.cos(k_phi) * np.sin(k_theta),
                                     np.sin(k_phi)], axis=1))
        h = np.cross(r, v_rel)
        expected_theta = np.degrees(np.arccos(h[:, 2] / np.linalg.norm(h, axis=1)))
        expected_phi = np.where(kick_phis > 0, mean_anomalies, mean_anomalies - 180.0)

        np.testing.assert_allclose(kick_info["theta_euler"], expected_theta, atol=1e-6)
        np.testing.assert_allclose(wrap_degrees(kick_info["phi_euler"] - expected_phi),
                                   0.0, atol=1e-6)

    def test_two_sn_systemic_velocity(self):
        """Test that the velocity changes from two SNe are combined correctly

        A He star + He star binary survives both SNe, with supplied kicks and
        mean anomalies for each. We independently follow the binary in the
        frame of the first SN (the frame the kick_info velocities are given
        in): the post-SN1 orbit has angular momentum direction h1 and
        periapsis direction e1 (for the separation r_1 - r_2), and COSMIC
        keeps that orbital plane and periapsis direction until the second SN
        (no precession). At the second SN, COSMIC's frame is built from
        r_2 - r_1, so its axes are x = -e1, z = h1, and the supplied kick
        angles are relative to those.

        The total change in systemic velocity must then match COSMIC's, which
        checks all three Euler angles used to rotate the second SN back into
        the first SN's frame (previously psi was off by a half-turn, and phi
        was random for circular pre-SN orbits).
        """
        G = 1.3271244e11  # km^3 / (Msun s^2), as in kick.f
        rsunkm = 6.96e5
        binary = InitialBinaryTable.InitialBinaries(
            m1=6.0, m2=5.5, porb=30.0, ecc=0.0, tphysf=3.0,
            kstar1=7, kstar2=7, metallicity=0.02
        )
        # [natal kick (km/s), phi (deg), theta (deg), mean anomaly (deg), seed]
        natal_kick_array = [[60.0, 40.0, 30.0, 100.0, 0.0],
                            [70.0, -25.0, 200.0, 250.0, 0.0]]
        bse = dict(BSEDict, kickflag=1, natal_kick_array=natal_kick_array)
        bpp, _, _, kick_info = Evolve.evolve(initialbinarytable=binary, BSEDict=bse,
                                             SSEDict=SSEDict, nproc=1, randomseed=42)
        bpp = bpp.reset_index(drop=True)
        kick_info = kick_info.set_index("star")
        self.assertTrue((kick_info["disrupted"] == 0).all())

        def kick_vector(vk, phi, theta, x, y, z):
            phi, theta = np.radians(phi), np.radians(theta)
            return vk * (np.cos(phi) * np.cos(theta) * x + np.cos(phi) * np.sin(theta) * y
                         + np.sin(phi) * z)

        def com_velocity_change(m_sn, m_rem, m_comp, w, kick):
            # velocities of each star relative to the pre-SN CoM, where w is
            # the relative velocity v_sn - v_comp, then the new CoM velocity
            m_tot = m_sn + m_comp
            v_sn, v_comp = m_comp / m_tot * w, -m_sn / m_tot * w
            return (m_rem * (v_sn + kick) + m_comp * v_comp) / (m_rem + m_comp)

        # ---- SN1 (star 1), in its own frame
        sn1 = bpp["evol_type"].tolist().index(15)
        m1, m2, a0, e0 = bpp.loc[sn1, ["mass_1", "mass_2", "sep", "ecc"]]
        m1_rem = bpp.loc[sn1 + 1, "mass_1"]
        s, w = kepler_state(a0 * rsunkm, e0, np.radians(natal_kick_array[0][3]), G * (m1 + m2))
        x, y, z = np.eye(3)
        k1 = kick_vector(*natal_kick_array[0][:3], x, y, z)
        V1 = com_velocity_change(m1, m1_rem, m2, w, k1)
        np.testing.assert_allclose(
            kick_info.loc[1, ["delta_vsysx_1", "delta_vsysy_1", "delta_vsysz_1"]].values,
            V1, atol=1e-6)

        # post-SN1 orbit of r_1 - r_2
        w1 = w + k1
        h1 = np.cross(s, w1)
        e1 = np.cross(w1, h1) / (G * (m1_rem + m2)) - s / np.linalg.norm(s)
        h1_hat, e1_hat = h1 / np.linalg.norm(h1), e1 / np.linalg.norm(e1)
        q1_hat = np.cross(h1_hat, e1_hat)

        # the stored Euler angles rotate (x, z) onto (e1, h1)
        R = euler_rotation(*kick_info.loc[1, ["theta_euler", "phi_euler", "psi_euler"]])
        np.testing.assert_allclose(R @ z, h1_hat, atol=1e-10)
        np.testing.assert_allclose(R @ x, e1_hat, atol=1e-10)

        # ---- SN2 (star 2), on the post-SN1 orbit with the pre-SN2 a, e, masses
        sn2 = bpp["evol_type"].tolist().index(16)
        m1b, m2b, a2, e2 = bpp.loc[sn2, ["mass_1", "mass_2", "sep", "ecc"]]
        m2_rem = bpp.loc[sn2 + 1, "mass_2"]
        s_in, w_in = kepler_state(a2 * rsunkm, e2, np.radians(natal_kick_array[1][3]),
                                  G * (m1b + m2b))
        # relative velocity of r_1 - r_2 in the SN1 frame
        w2 = w_in[0] * e1_hat + w_in[1] * q1_hat
        # COSMIC's SN2 frame (separation r_2 - r_1): x = -e1, z = h1
        x2, z2 = -e1_hat, h1_hat
        k2 = kick_vector(*natal_kick_array[1][:3], x2, np.cross(z2, x2), z2)
        V2 = com_velocity_change(m2b, m2_rem, m1b, -w2, k2)

        V_total = V1 + V2
        cosmic_total = sum(kick_info.loc[star, ["delta_vsysx_1", "delta_vsysy_1",
                                                "delta_vsysz_1"]].values for star in [1, 2])
        np.testing.assert_allclose(cosmic_total, V_total, atol=1e-6)
        self.assertAlmostEqual(kick_info.loc[2, "vsys_1_total"], np.linalg.norm(V_total),
                               delta=1e-6)

    def test_two_sn_disrupted_runaway_velocities(self):
        """Test the stars' total velocities when the second SN disrupts the binary

        The binary survives the first SN (with a small kick) and is disrupted
        by the second (with a large one). Each star's final velocity, relative
        to the original centre of mass and in the first SN's frame, is then
        the systemic velocity from the first SN plus the change in centre of
        mass velocity from the second, plus its share of the relative velocity
        at infinity of the hyperbolic post-SN orbit, which we find
        independently in its perifocal frame.

        Previously vsys_2_total (the total speed of star 2, which explodes
        second) left out the systemic velocity from the first SN.
        """
        G = 1.3271244e11  # km^3 / (Msun s^2), as in kick.f
        rsunkm = 6.96e5
        binary = InitialBinaryTable.InitialBinaries(
            m1=6.0, m2=5.5, porb=30.0, ecc=0.0, tphysf=3.0,
            kstar1=7, kstar2=7, metallicity=0.02
        )
        # [natal kick (km/s), phi (deg), theta (deg), mean anomaly (deg), seed]
        natal_kick_array = [[60.0, 40.0, 30.0, 100.0, 0.0],
                            [700.0, -25.0, 200.0, 250.0, 0.0]]
        bse = dict(BSEDict, kickflag=1, natal_kick_array=natal_kick_array)
        bpp, _, _, kick_info = Evolve.evolve(initialbinarytable=binary, BSEDict=bse,
                                             SSEDict=SSEDict, nproc=1, randomseed=42)
        bpp = bpp.reset_index(drop=True)
        kick_info = kick_info.set_index("star")
        self.assertEqual(kick_info.loc[1, "disrupted"], 0)
        self.assertEqual(kick_info.loc[2, "disrupted"], 1)

        def kick_vector(vk, phi, theta, x, y, z):
            phi, theta = np.radians(phi), np.radians(theta)
            return vk * (np.cos(phi) * np.cos(theta) * x + np.cos(phi) * np.sin(theta) * y
                         + np.sin(phi) * z)

        def com_velocity_change(m_sn, m_rem, m_comp, w, kick):
            m_tot = m_sn + m_comp
            v_sn, v_comp = m_comp / m_tot * w, -m_sn / m_tot * w
            return (m_rem * (v_sn + kick) + m_comp * v_comp) / (m_rem + m_comp)

        # ---- SN1 (star 1), which leaves the binary bound
        sn1 = bpp["evol_type"].tolist().index(15)
        m1, m2, a0, e0 = bpp.loc[sn1, ["mass_1", "mass_2", "sep", "ecc"]]
        m1_rem = bpp.loc[sn1 + 1, "mass_1"]
        s, w = kepler_state(a0 * rsunkm, e0, np.radians(natal_kick_array[0][3]), G * (m1 + m2))
        x, y, z = np.eye(3)
        k1 = kick_vector(*natal_kick_array[0][:3], x, y, z)
        V1 = com_velocity_change(m1, m1_rem, m2, w, k1)
        w1 = w + k1
        h1 = np.cross(s, w1)
        e1 = np.cross(w1, h1) / (G * (m1_rem + m2)) - s / np.linalg.norm(s)
        h1_hat, e1_hat = h1 / np.linalg.norm(h1), e1 / np.linalg.norm(e1)
        q1_hat = np.cross(h1_hat, e1_hat)

        # ---- SN2 (star 2), which disrupts it
        sn2 = bpp["evol_type"].tolist().index(16)
        m1b, m2b, a2, e2 = bpp.loc[sn2, ["mass_1", "mass_2", "sep", "ecc"]]
        m2_rem = bpp.loc[sn2 + 1, "mass_2"]
        s_in, w_in = kepler_state(a2 * rsunkm, e2, np.radians(natal_kick_array[1][3]),
                                  G * (m1b + m2b))
        # separation r_1 - r_2 and its rate of change, in the SN1 frame
        s2 = s_in[0] * e1_hat + s_in[1] * q1_hat
        w2 = w_in[0] * e1_hat + w_in[1] * q1_hat
        x2, z2 = -e1_hat, h1_hat
        k2 = kick_vector(*natal_kick_array[1][:3], x2, np.cross(z2, x2), z2)
        V2 = com_velocity_change(m2b, m2_rem, m1b, -w2, k2)

        # hyperbolic orbit of r_2 - r_1 after SN2
        gm = G * (m2_rem + m1b)
        r_rel, v_rel = -s2, -w2 + k2
        h = np.cross(r_rel, v_rel)
        e_vec = np.cross(v_rel, h) / gm - r_rel / np.linalg.norm(r_rel)
        ecc = np.linalg.norm(e_vec)
        self.assertGreater(ecc, 1.0)
        # speed at infinity from energy, direction from the outgoing asymptote
        # (true anomaly arccos(-1/e)), where the perifocal velocity is
        # proportional to (-sin nu, e + cos nu)
        v_inf = np.sqrt(v_rel @ v_rel - 2 * gm / np.linalg.norm(r_rel))
        nu_inf = np.arccos(-1 / ecc)
        p_hat = e_vec / ecc
        q_hat = np.cross(h / np.linalg.norm(h), p_hat)
        u = -np.sin(nu_inf) * p_hat + (ecc + np.cos(nu_inf)) * q_hat
        v_inf_vec = v_inf * u / np.linalg.norm(u)

        v_star2 = V1 + V2 + m1b / (m2_rem + m1b) * v_inf_vec
        v_star1 = V1 + V2 - m2_rem / (m2_rem + m1b) * v_inf_vec

        sn1_vsys = kick_info.loc[1, ["delta_vsysx_1", "delta_vsysy_1", "delta_vsysz_1"]].values
        np.testing.assert_allclose(sn1_vsys, V1, atol=1e-6)
        np.testing.assert_allclose(
            sn1_vsys + kick_info.loc[2, ["delta_vsysx_1", "delta_vsysy_1",
                                         "delta_vsysz_1"]].values, v_star1, atol=1e-6)
        np.testing.assert_allclose(
            sn1_vsys + kick_info.loc[2, ["delta_vsysx_2", "delta_vsysy_2",
                                         "delta_vsysz_2"]].values, v_star2, atol=1e-6)
        self.assertAlmostEqual(kick_info.loc[2, "vsys_1_total"], np.linalg.norm(v_star1),
                               delta=1e-6)
        self.assertAlmostEqual(kick_info.loc[2, "vsys_2_total"], np.linalg.norm(v_star2),
                               delta=1e-6)

    def test_psi_euler_periapsis_on_line_of_nodes(self):
        """Test psi_euler when the new periapsis lies on the line of nodes

        The line of nodes always passes through the stars' separation at the
        SN, so if the kick leaves no radial velocity the star is at the new
        periapsis (or apoapsis), on the line of nodes, and psi_euler is 0 (or
        180 deg). In a circular orbit that happens for a kick at
        theta = M + 90 deg. The cosine of psi can then round to just above 1,
        and arccos used to return NaN (about 1 time in 5).
        """
        binary = InitialBinaryTable.InitialBinaries(
            m1=6.0, m2=10.0, porb=5.0, ecc=0.0, tphysf=2.0,
            kstar1=7, kstar2=1, metallicity=0.02
        )
        mean_anomalies = np.arange(0.0, 360.0, 10.0)
        n = len(mean_anomalies)
        bpp, kick_info = evolve_with_first_sn_kicks(
            binary, natal_kick=np.full(n, 50.0), phi=np.full(n, 40.0),
            theta=(mean_anomalies + 90.0) % 360.0, mean_anomaly=mean_anomalies, tphysf=1.2
        )
        self.assertTrue((bpp[bpp["evol_type"] == 15]["ecc"] == 0.0).all())
        self.assertTrue((kick_info["disrupted"] == 0).all())
        psi = kick_info["psi_euler"].values
        self.assertTrue(np.isfinite(psi).all())
        np.testing.assert_allclose(np.sin(np.radians(psi)), 0.0, atol=1e-6)

    def test_supplied_kick_magnitude_random_stream(self):
        """Test that supplying a kick speed doesn't change the later random draws

        When a natal kick speed is supplied, COSMIC calls the random number
        generator as many times as drawing the speed would have, so that the
        kick direction (and anything else drawn later) is the same as if the
        speed had been drawn. Previously this was only right for the Maxwellian
        kicks, and the truncated distributions (kickflag 6 and 9) used a
        varying number of draws.

        Here single stars that form NSs and BHs are evolved once drawing their
        kick speeds, then again (with the same seeds) supplying those speeds.
        """
        n = 40
        stars = pd.concat([InitialBinaryTable.InitialBinaries(
            m1=m1, m2=0.0, porb=0.0, ecc=-1.0, tphysf=100.0,
            kstar1=1, kstar2=15, metallicity=0.02
        ) for m1 in [20.005, 30.0]])
        stars = stars.loc[stars.index.repeat(n)].reset_index(drop=True)
        stars["randomseed"] = -(np.arange(len(stars)) + 1234)
        bse = {k: v for k, v in BSEDict.items()
               if k != "natal_kick_array" and k not in FLATTENED_NATAL_KICK_COLUMNS}
        bse.update(ecsn=0.0, ecsn_mlow=0.0, bhflag=1)
        for kickflag in [1, 5, 6, 9]:
            with self.subTest(kickflag=kickflag):
                bse["kickflag"] = kickflag
                tables = []
                for supplied in [None, "previous"]:
                    table = stars.copy()
                    for col in ["natal_kick", "phi", "theta", "mean_anomaly"]:
                        table[col + "_1"] = -100.0
                        table[col + "_2"] = -100.0
                    table["randomseed_1"] = 0.0
                    table["randomseed_2"] = 0.0
                    if supplied is not None:
                        table["natal_kick_1"] = tables[0]["natal_kick"].values
                    _, _, _, kick_info = Evolve.evolve(initialbinarytable=table,
                                                       BSEDict=bse, SSEDict=SSEDict,
                                                       nproc=1)
                    tables.append(kick_info[kick_info["star"] == 1]
                                  .sort_index()[["natal_kick", "phi", "theta"]])
                drawn, supplied = tables
                self.assertEqual(len(drawn), len(stars))
                # (BH kicks can be zero after fallback)
                self.assertTrue((drawn["natal_kick"] > 0).any())
                pd.testing.assert_frame_equal(drawn, supplied)
