"""Unit tests for the STROOPWAFEL adaptive importance sampling module.

Covers the composable distribution/transform design
(`cosmic.sample.stroopwafel.distributions`), the vectorized
`ParameterSpace`, physical rejection, the Gaussian mixture model, and the
`COSMICStroopOutput` result type.
"""

__author__ = 'Tom Wagg <tom.wagg@gmail.com>'

import unittest
import numpy as np
import pandas as pd
from scipy.stats import norm, kstest
from scipy.integrate import trapezoid

from cosmic.sample.stroopwafel import ParameterSpace, Parameter
from cosmic.sample.stroopwafel.distributions import (
    DISTRIBUTIONS, register, get_distribution,
    Distribution, Uniform, PowerLaw, TruncatedNormal,
    Identity, Log10, Ln, Sin, CosShift,
)
from cosmic.sample.stroopwafel.mixture_model import GaussianMixture
from cosmic.sample.stroopwafel.rejection import get_zams_radius, default_reject
from cosmic.sample.stroopwafel.constants import (
    ALPHA_IMF, SANA_G, SANA_ECC, NATAL_KICK_LOG_MU, NATAL_KICK_LOG_SIGMA,
)
from cosmic.output import COSMICStroopOutput

HALF_PI = np.pi / 2
KS_PVALUE_MIN = 0.01   # samplers must not be rejected against their own CDF


# --------------------------------------------------------------------------
# Theoretical CDFs used for goodness-of-fit tests
# --------------------------------------------------------------------------
def _uniform_cdf(x, lo, hi):
    return (x - lo) / (hi - lo)


def _powerlaw_cdf(x, alpha, lo, hi):
    a = alpha + 1
    return (np.power(x, a) - lo**a) / (hi**a - lo**a)


def _truncnorm_cdf(x, mu, scale, lo, hi):
    lo_c, hi_c = norm.cdf(lo, mu, scale), norm.cdf(hi, mu, scale)
    return (norm.cdf(x, mu, scale) - lo_c) / (hi_c - lo_c)


def canonical_space():
    """The reference 6-D parameter space exercising every built-in dist."""
    return ParameterSpace([
        Parameter('mass_1', 5.0, 150.0, dist='kroupa'),
        Parameter('q', 0.0, 1.0, dist='uniform'),
        Parameter('porb', 10**(0.15), 10**(5.5), dist='sana'),
        Parameter('ecc', 1e-9, 0.99999999, dist='sana_ecc'),
        Parameter('metallicity', 1e-4, 0.03, dist='flat_in_log'),
        Parameter('natal_kick_1', 0.1, 5000.0, dist='log_normal'),
    ])


class Triangular(Distribution):
    """A user-defined distribution, used to test extensibility."""

    def sample(self, n, lo, hi, rng=None):
        rng = rng or np.random.default_rng()
        return rng.triangular(lo, 0.5 * (lo + hi), hi, n)

    def pdf(self, values, lo, hi):
        mode = 0.5 * (lo + hi)
        height = 2.0 / (hi - lo)
        return np.where(values < mode,
                        height * (values - lo) / (mode - lo),
                        height * (hi - values) / (hi - mode))


# ==========================================================================
# Transforms
# ==========================================================================
class TestTransforms(unittest.TestCase):

    def test_identity(self):
        t = Identity()
        x = np.array([-3.0, 0.0, 7.5])
        np.testing.assert_array_equal(t.to_sampling(x), x)
        np.testing.assert_array_equal(t.to_physical(x), x)
        self.assertEqual(t.bounds(2.0, 9.0), (2.0, 9.0))

    def test_log10(self):
        t = Log10()
        x = np.geomspace(1e-4, 1e3, 50)
        np.testing.assert_allclose(t.to_physical(t.to_sampling(x)), x, rtol=1e-12)
        lo, hi = t.bounds(1e-4, 1e3)
        self.assertAlmostEqual(lo, -4.0)
        self.assertAlmostEqual(hi, 3.0)

    def test_ln(self):
        t = Ln()
        x = np.geomspace(0.1, 5000.0, 50)
        np.testing.assert_allclose(t.to_physical(t.to_sampling(x)), x, rtol=1e-12)
        lo, hi = t.bounds(0.1, 5000.0)
        self.assertAlmostEqual(lo, np.log(0.1))
        self.assertAlmostEqual(hi, np.log(5000.0))

    def test_sin(self):
        t = Sin()
        x = np.linspace(-HALF_PI, HALF_PI, 50)
        np.testing.assert_allclose(t.to_physical(t.to_sampling(x)), x, atol=1e-12)
        self.assertEqual(t.bounds(-HALF_PI, HALF_PI), (-1.0, 1.0))

    def test_cosshift_bounds_are_sorted(self):
        # CosShift is monotonically decreasing, so bounds() must swap ends.
        t = CosShift()
        x = np.linspace(-HALF_PI, HALF_PI, 50)
        np.testing.assert_allclose(t.to_physical(t.to_sampling(x)), x, atol=1e-12)
        lo, hi = t.bounds(-HALF_PI, HALF_PI)
        self.assertLess(lo, hi)
        np.testing.assert_allclose([lo, hi], [-1.0, 1.0], atol=1e-12)


# ==========================================================================
# Distributions
# ==========================================================================
class TestDistributions(unittest.TestCase):

    def test_registry_contents(self):
        self.assertEqual(
            set(DISTRIBUTIONS),
            {'uniform', 'flat_in_log', 'uniform_in_sine', 'uniform_in_cosine',
             'kroupa', 'sana', 'sana_ecc', 'log_normal'},
        )

    def test_builtin_compositions(self):
        self.assertIsInstance(DISTRIBUTIONS['kroupa'], PowerLaw)
        self.assertIsInstance(DISTRIBUTIONS['kroupa'].transform, Identity)
        self.assertIsInstance(DISTRIBUTIONS['sana'], PowerLaw)
        self.assertIsInstance(DISTRIBUTIONS['sana'].transform, Log10)
        self.assertIsInstance(DISTRIBUTIONS['flat_in_log'], Uniform)
        self.assertIsInstance(DISTRIBUTIONS['flat_in_log'].transform, Log10)
        self.assertIsInstance(DISTRIBUTIONS['log_normal'], TruncatedNormal)
        self.assertIsInstance(DISTRIBUTIONS['log_normal'].transform, Ln)
        self.assertEqual(DISTRIBUTIONS['kroupa'].alpha, ALPHA_IMF)
        self.assertEqual(DISTRIBUTIONS['sana'].alpha, SANA_G)
        self.assertEqual(DISTRIBUTIONS['sana_ecc'].alpha, SANA_ECC)

    def test_uniform_pdf_constant(self):
        d = Uniform()
        vals = np.linspace(2.0, 8.0, 100)
        np.testing.assert_allclose(d.pdf(vals, 2.0, 8.0), 1.0 / 6.0)

    def test_powerlaw_pdf_matches_closed_form(self):
        for alpha, lo, hi in [(ALPHA_IMF, 5.0, 150.0),
                              (SANA_G, 0.15, 5.5),
                              (SANA_ECC, 1e-9, 0.999)]:
            vals = np.linspace(lo * 1.01, hi * 0.99, 100)
            got = PowerLaw(alpha).pdf(vals, lo, hi)
            a = alpha + 1
            expected = (a / (hi**a - lo**a)) * vals**alpha
            np.testing.assert_allclose(got, expected, rtol=1e-12)

    def test_powerlaw_pdf_normalised(self):
        d = PowerLaw(ALPHA_IMF)
        lo, hi = 5.0, 150.0
        x = np.linspace(lo, hi, 200000)
        self.assertAlmostEqual(trapezoid(d.pdf(x, lo, hi), x), 1.0, places=4)

    def test_truncnorm_pdf_normalised(self):
        d = TruncatedNormal(NATAL_KICK_LOG_MU, NATAL_KICK_LOG_SIGMA)
        lo, hi = 2.0, 9.0
        x = np.linspace(lo, hi, 200000)
        self.assertAlmostEqual(trapezoid(d.pdf(x, lo, hi), x), 1.0, places=4)

    def test_uniform_sampler_ks(self):
        rng = np.random.default_rng(0)
        s = Uniform().sample(20000, 2.0, 8.0, rng=rng)
        self.assertTrue(np.all((s >= 2.0) & (s <= 8.0)))
        p = kstest(s, lambda v: _uniform_cdf(v, 2.0, 8.0)).pvalue
        self.assertGreater(p, KS_PVALUE_MIN)

    def test_powerlaw_sampler_ks(self):
        for alpha, lo, hi, seed in [(ALPHA_IMF, 5.0, 150.0, 1),
                                    (SANA_G, 0.15, 5.5, 2),
                                    (SANA_ECC, 1e-9, 0.999, 3)]:
            rng = np.random.default_rng(seed)
            s = PowerLaw(alpha).sample(20000, lo, hi, rng=rng)
            self.assertTrue(np.all((s >= lo) & (s <= hi)))
            p = kstest(s, lambda v: _powerlaw_cdf(v, alpha, lo, hi)).pvalue
            self.assertGreater(p, KS_PVALUE_MIN, msg=f"alpha={alpha}")

    def test_truncnorm_sampler_ks(self):
        rng = np.random.default_rng(4)
        mu, scale, lo, hi = NATAL_KICK_LOG_MU, NATAL_KICK_LOG_SIGMA, 2.0, 9.0
        s = TruncatedNormal(mu, scale).sample(20000, lo, hi, rng=rng)
        self.assertTrue(np.all((s >= lo) & (s <= hi)))
        p = kstest(s, lambda v: _truncnorm_cdf(v, mu, scale, lo, hi)).pvalue
        self.assertGreater(p, KS_PVALUE_MIN)

    def test_default_sigma_is_avg_density_over_pdf(self):
        d = Uniform()
        vals = np.linspace(2.0, 8.0, 50)
        ad = 0.05
        np.testing.assert_allclose(
            d.sigma(vals, 2.0, 8.0, ad), ad / d.pdf(vals, 2.0, 8.0), rtol=1e-12,
        )

    def test_powerlaw_sigma_positive_and_raises_on_nonpositive_lo(self):
        d = PowerLaw(ALPHA_IMF)
        vals = np.linspace(6.0, 140.0, 100)
        sig = d.sigma(vals, 5.0, 150.0, 0.05)
        self.assertTrue(np.all(sig > 0) and np.all(np.isfinite(sig)))
        with self.assertRaises(ValueError):
            d.sigma(vals, 0.0, 150.0, 0.05)


# ==========================================================================
# Registry / extensibility
# ==========================================================================
class TestRegistry(unittest.TestCase):

    def test_get_distribution_passthrough(self):
        d = PowerLaw(-1.7)
        self.assertIs(get_distribution(d), d)

    def test_get_distribution_unknown_raises(self):
        with self.assertRaises(KeyError):
            get_distribution('does_not_exist')

    def test_register_and_lookup(self):
        self.addCleanup(lambda: DISTRIBUTIONS.pop('test_steep_imf', None))
        register('test_steep_imf', PowerLaw(-2.7))
        p = Parameter('m', 5.0, 100.0, dist='test_steep_imf')
        self.assertEqual(p.distribution.alpha, -2.7)

    def test_register_rejects_non_instance(self):
        with self.assertRaises(TypeError):
            register('bad', PowerLaw)   # a class, not an instance

    def test_custom_distribution_instance(self):
        p = Parameter('x', 1.0, 10.0, dist=PowerLaw(-1.7))
        s = p.distribution.sample(500, p.lo, p.hi, rng=np.random.default_rng(0))
        self.assertTrue(np.all((s >= 1.0) & (s <= 10.0)))

    def test_custom_subclass_end_to_end(self):
        # A user-defined Distribution composed with a transform, run through
        # the full ParameterSpace pipeline.
        ps = ParameterSpace([Parameter('t', 1.0, 100.0, dist=Triangular(transform=Log10()))])
        samples, mask = ps.sample(5000, rng=np.random.default_rng(1))
        phys = ps.to_physical(samples[mask])
        self.assertTrue(np.all((phys >= 1.0) & (phys <= 100.0)))
        self.assertTrue(np.all(ps.compute_prior(samples[mask]) > 0))
        sig = ps.compute_sigma(samples[mask][:100], 0.05)
        self.assertTrue(np.all(sig > 0) and np.all(np.isfinite(sig)))


# ==========================================================================
# ParameterSpace
# ==========================================================================
class TestParameterSpace(unittest.TestCase):

    def setUp(self):
        self.ps = canonical_space()

    def test_names_sorted_and_ndim(self):
        self.assertEqual(self.ps.names, sorted(self.ps.names))
        self.assertEqual(self.ps.ndim, 6)

    def test_idx(self):
        for i, name in enumerate(self.ps.names):
            self.assertEqual(self.ps.idx(name), i)

    def test_sample_shapes_and_mask(self):
        samples, mask = self.ps.sample(1000, rng=np.random.default_rng(42))
        self.assertEqual(samples.shape, (1000, 6))
        self.assertEqual(mask.shape, (1000,))
        self.assertEqual(mask.dtype, bool)

    def test_roundtrip_transform(self):
        samples, _ = self.ps.sample(1000, rng=np.random.default_rng(42))
        back = self.ps.to_sampling(self.ps.to_physical(samples))
        np.testing.assert_allclose(back, samples, atol=1e-12)

    def test_prior_positive_finite(self):
        samples, mask = self.ps.sample(1000, rng=np.random.default_rng(42))
        priors = self.ps.compute_prior(samples[mask])
        self.assertTrue(np.all(priors > 0) and np.all(np.isfinite(priors)))

    def test_in_bounds_matches_sample_mask(self):
        samples, mask = self.ps.sample(1000, rng=np.random.default_rng(42))
        np.testing.assert_array_equal(self.ps.in_bounds(samples), mask)

    def test_sana_bounds_linear_to_log10(self):
        # sana bounds are given in linear (physical) days; the Log10 transform
        # maps them into log10-period sampling space.
        porb = next(p for p in self.ps.params if p.name == 'porb')
        self.assertAlmostEqual(porb.lo, 0.15)
        self.assertAlmostEqual(porb.hi, 5.5)

    def test_compute_sigma_positive_finite(self):
        samples, mask = self.ps.sample(5000, rng=np.random.default_rng(7))
        hits = samples[mask][:200]
        ad = 1.0 / np.power(5000, 1.0 / self.ps.ndim)
        sig = self.ps.compute_sigma(hits, ad)
        self.assertEqual(sig.shape, hits.shape)
        self.assertTrue(np.all(sig > 0) and np.all(np.isfinite(sig)))

    def test_compute_sigma_reports_parameter_name_on_error(self):
        # A sana period with a 1-day lower bound maps to log10(1)=0, which the
        # power-law sigma cannot handle; the error should name the parameter.
        ps = ParameterSpace([Parameter('porb', 1.0, 1000.0, dist='sana')])
        samples, mask = ps.sample(100, rng=np.random.default_rng(0))
        with self.assertRaisesRegex(ValueError, 'porb'):
            ps.compute_sigma(samples[mask][:10], 0.05)


# ==========================================================================
# Rejection
# ==========================================================================
class TestRejection(unittest.TestCase):

    PARAM_NAMES = ['ecc', 'mass_1', 'metallicity', 'porb', 'q']

    def test_get_zams_radius_positive_finite(self):
        masses = np.array([1.0, 10.0, 30.0, 100.0])
        mets = np.array([0.02, 0.02, 0.001, 0.014])
        r = get_zams_radius(masses, mets)
        self.assertEqual(r.shape, (4,))
        self.assertTrue(np.all(r > 0) and np.all(np.isfinite(r)))

    def test_wide_binary_not_rejected(self):
        samples = np.array([[0.1, 20.0, 0.014, 100.0, 0.5]])
        derived = {
            'mass_2': np.array([10.0]),
            'metallicity_1': np.array([0.014]),
            'metallicity_2': np.array([0.014]),
            'separation': np.array([50.0]),
        }
        rejected = default_reject(samples, derived, self.PARAM_NAMES)
        self.assertFalse(rejected[0])

    def test_low_mass_secondary_rejected(self):
        samples = np.array([[0.1, 20.0, 0.014, 100.0, 0.001]])
        derived = {
            'mass_2': np.array([0.02]),   # below the 0.08 Msun minimum
            'metallicity_1': np.array([0.014]),
            'metallicity_2': np.array([0.014]),
            'separation': np.array([50.0]),
        }
        rejected = default_reject(samples, derived, self.PARAM_NAMES)
        self.assertTrue(rejected[0])


# ==========================================================================
# Gaussian mixture model
# ==========================================================================
class TestGaussianMixture(unittest.TestCase):

    def setUp(self):
        self.ps = canonical_space()
        rng = np.random.default_rng(42)
        samples, mask = self.ps.sample(20000, rng=rng)
        self.valid = samples[mask]
        self.hits = self.valid[:10]
        self.ad = 1.0 / np.power(20000, 1.0 / self.ps.ndim)

    def test_from_hits_shapes(self):
        gm = GaussianMixture.from_hits(self.hits, self.ps, self.ad)
        self.assertEqual(gm.n_components, 10)
        self.assertEqual(gm.means.shape, (10, 6))
        self.assertEqual(gm.covariances.shape, (10, 6, 6))
        self.assertAlmostEqual(float(np.sum(gm.alphas)), 1.0)

    def test_pdf_nonnegative_finite(self):
        gm = GaussianMixture.from_hits(self.hits, self.ps, self.ad)
        pdf = gm.pdf(self.valid[:200])
        self.assertTrue(np.all(pdf >= 0) and np.all(np.isfinite(pdf)))

    def test_sample_shapes_and_idx_range(self):
        gm = GaussianMixture.from_hits(self.hits, self.ps, self.ad)
        msamp, mmask, midx = gm.sample(2000, self.ps, rng=np.random.default_rng(1))
        self.assertEqual(msamp.shape[1], 6)
        self.assertEqual(len(mmask), len(msamp))
        self.assertEqual(len(midx), len(msamp))
        self.assertTrue(np.all((midx >= 0) & (midx < 10)))


# ==========================================================================
# COSMICStroopOutput
# ==========================================================================
class TestStroopOutput(unittest.TestCase):

    def test_hit_rate(self):
        N = 100
        bin_nums = np.arange(N)
        is_hit = np.zeros(N, dtype=bool)
        is_hit[:10] = True
        frame = pd.DataFrame({'bin_num': bin_nums})
        result = COSMICStroopOutput(
            bpp=frame, bcm=frame, initC=frame, kick_info=frame,
            samples=np.zeros((N, 1)), param_names=['mass_1'],
            weights=np.ones(N), is_hit=is_hit,
            generation=np.zeros(N, dtype=int),
            gaussian_idx=np.full(N, -1, dtype=int),
            num_explored=N, num_hits=10, fraction_explored=1.0,
        )
        self.assertAlmostEqual(result.hit_rate, 0.1)


if __name__ == '__main__':
    unittest.main()
