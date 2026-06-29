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

from cosmic.sample.stroopwafel import ParameterSpace, Parameter, AdaptiveSampler
from cosmic.sample.stroopwafel.distributions import (
    DISTRIBUTIONS, register, get_distribution,
    Distribution, Uniform, PowerLaw, BrokenPowerLaw, TruncatedNormal,
    Identity, Log10, Ln, Sin, CosShift,
)
from cosmic.sample.stroopwafel.mixture_model import GaussianMixture
from cosmic.sample.stroopwafel.rejection import default_reject
from scipy.integrate import cumulative_trapezoid
from cosmic.output import COSMICStroopOutput

# Built-in distribution parameters, kept in sync with the DISTRIBUTIONS registry.
KROUPA_BREAK, KROUPA_ALPHA_LOW, KROUPA_ALPHA_HIGH = 0.5, -1.3, -2.3
SANA_G, SANA_ECC = -0.55, -0.45
DISBERG_MU, DISBERG_SIGMA = 5.67, 0.59

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
        Parameter('natal_kick_1', 0.1, 5000.0, dist='disberg'),
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
             'kroupa', 'sana', 'sana_ecc', 'disberg'},
        )

    def test_builtin_compositions(self):
        self.assertIsInstance(DISTRIBUTIONS['kroupa'], BrokenPowerLaw)
        self.assertIsInstance(DISTRIBUTIONS['kroupa'].transform, Identity)
        np.testing.assert_array_equal(DISTRIBUTIONS['kroupa'].breaks, [KROUPA_BREAK])
        np.testing.assert_array_equal(DISTRIBUTIONS['kroupa'].alphas,
                                      [KROUPA_ALPHA_LOW, KROUPA_ALPHA_HIGH])
        self.assertIsInstance(DISTRIBUTIONS['sana'], PowerLaw)
        self.assertIsInstance(DISTRIBUTIONS['sana'].transform, Log10)
        self.assertIsInstance(DISTRIBUTIONS['flat_in_log'], Uniform)
        self.assertIsInstance(DISTRIBUTIONS['flat_in_log'].transform, Log10)
        self.assertIsInstance(DISTRIBUTIONS['disberg'], TruncatedNormal)
        self.assertIsInstance(DISTRIBUTIONS['disberg'].transform, Ln)
        self.assertEqual(DISTRIBUTIONS['sana'].alpha, SANA_G)
        self.assertEqual(DISTRIBUTIONS['sana_ecc'].alpha, SANA_ECC)

    def test_uniform_pdf_constant(self):
        d = Uniform()
        vals = np.linspace(2.0, 8.0, 100)
        np.testing.assert_allclose(d.pdf(vals, 2.0, 8.0), 1.0 / 6.0)

    def test_powerlaw_pdf_matches_closed_form(self):
        for alpha, lo, hi in [(KROUPA_ALPHA_HIGH, 5.0, 150.0),
                              (SANA_G, 0.15, 5.5),
                              (SANA_ECC, 1e-9, 0.999)]:
            vals = np.linspace(lo * 1.01, hi * 0.99, 100)
            got = PowerLaw(alpha).pdf(vals, lo, hi)
            a = alpha + 1
            expected = (a / (hi**a - lo**a)) * vals**alpha
            np.testing.assert_allclose(got, expected, rtol=1e-12)

    def test_powerlaw_pdf_normalised(self):
        d = PowerLaw(KROUPA_ALPHA_HIGH)
        lo, hi = 5.0, 150.0
        x = np.linspace(lo, hi, 200000)
        self.assertAlmostEqual(trapezoid(d.pdf(x, lo, hi), x), 1.0, places=4)

    def test_truncnorm_pdf_normalised(self):
        d = TruncatedNormal(DISBERG_MU, DISBERG_SIGMA)
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
        for alpha, lo, hi, seed in [(KROUPA_ALPHA_HIGH, 5.0, 150.0, 1),
                                    (SANA_G, 0.15, 5.5, 2),
                                    (SANA_ECC, 1e-9, 0.999, 3)]:
            rng = np.random.default_rng(seed)
            s = PowerLaw(alpha).sample(20000, lo, hi, rng=rng)
            self.assertTrue(np.all((s >= lo) & (s <= hi)))
            p = kstest(s, lambda v: _powerlaw_cdf(v, alpha, lo, hi)).pvalue
            self.assertGreater(p, KS_PVALUE_MIN, msg=f"alpha={alpha}")

    def test_truncnorm_sampler_ks(self):
        rng = np.random.default_rng(4)
        mu, scale, lo, hi = DISBERG_MU, DISBERG_SIGMA, 2.0, 9.0
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
        d = PowerLaw(KROUPA_ALPHA_HIGH)
        vals = np.linspace(6.0, 140.0, 100)
        sig = d.sigma(vals, 5.0, 150.0, 0.05)
        self.assertTrue(np.all(sig > 0) and np.all(np.isfinite(sig)))
        with self.assertRaises(ValueError):
            d.sigma(vals, 0.0, 150.0, 0.05)


# ==========================================================================
# Broken power law (Kroupa IMF)
# ==========================================================================
class TestBrokenPowerLaw(unittest.TestCase):

    KROUPA = BrokenPowerLaw(breaks=[KROUPA_BREAK], alphas=[KROUPA_ALPHA_LOW, KROUPA_ALPHA_HIGH])

    def test_input_validation(self):
        with self.assertRaises(ValueError):
            BrokenPowerLaw(breaks=[0.5], alphas=[-2.3])              # wrong length
        with self.assertRaises(ValueError):
            BrokenPowerLaw(breaks=[0.5, 0.3], alphas=[-1, -2, -3])   # not increasing
        with self.assertRaises(ValueError):
            BrokenPowerLaw(breaks=[0.5], alphas=[-1.0, -2.3])        # exponent -1

    def test_pdf_continuous_at_break(self):
        lo, hi = 0.08, 100.0
        below = self.KROUPA.pdf(np.array([KROUPA_BREAK - 1e-6]), lo, hi)[0]
        above = self.KROUPA.pdf(np.array([KROUPA_BREAK + 1e-6]), lo, hi)[0]
        self.assertAlmostEqual(below, above, places=4)

    def test_pdf_normalised_across_break(self):
        lo, hi = 0.08, 100.0
        x = np.geomspace(lo, hi, 200000)
        self.assertAlmostEqual(trapezoid(self.KROUPA.pdf(x, lo, hi), x), 1.0, places=3)

    def test_slope_changes_at_break(self):
        # Local log-log slope should match the segment exponents.
        lo, hi = 0.08, 100.0
        def local_slope(x0):
            x = np.array([x0 * 0.999, x0 * 1.001])
            p = self.KROUPA.pdf(x, lo, hi)
            return np.diff(np.log(p))[0] / np.diff(np.log(x))[0]
        self.assertAlmostEqual(local_slope(0.2), KROUPA_ALPHA_LOW, places=2)
        self.assertAlmostEqual(local_slope(5.0), KROUPA_ALPHA_HIGH, places=2)

    def test_reduces_to_powerlaw_above_break(self):
        # With no break in range the example regime (m1 > 5) is unchanged.
        lo, hi = 5.0, 150.0
        pl = PowerLaw(KROUPA_ALPHA_HIGH)
        vals = np.linspace(lo * 1.01, hi * 0.99, 200)
        np.testing.assert_allclose(self.KROUPA.pdf(vals, lo, hi),
                                   pl.pdf(vals, lo, hi), rtol=1e-10)
        a = self.KROUPA.sample(2000, lo, hi, rng=np.random.default_rng(0))
        b = pl.sample(2000, lo, hi, rng=np.random.default_rng(0))
        np.testing.assert_allclose(a, b, rtol=1e-10)

    def test_sampler_follows_pdf(self):
        # KS test against the (independently integrated) pdf, across the break.
        lo, hi = 0.08, 50.0
        s = self.KROUPA.sample(40000, lo, hi, rng=np.random.default_rng(7))
        self.assertTrue(np.all((s >= lo) & (s <= hi)))
        grid = np.geomspace(lo, hi, 40000)
        cdf_vals = cumulative_trapezoid(self.KROUPA.pdf(grid, lo, hi), grid, initial=0)
        cdf_vals /= cdf_vals[-1]
        pvalue = kstest(s, lambda v: np.interp(v, grid, cdf_vals)).pvalue
        self.assertGreater(pvalue, KS_PVALUE_MIN)

    def test_sigma_positive_finite_across_break(self):
        lo, hi = 0.08, 100.0
        vals = np.array([0.1, 0.3, 0.5, 1.0, 10.0, 80.0])
        sig = self.KROUPA.sigma(vals, lo, hi, 0.02)
        self.assertTrue(np.all(sig > 0) and np.all(np.isfinite(sig)))


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

    def test_wide_binary_not_rejected(self):
        binary_params = {
            'mass_1': np.array([20.0]),
            'mass_2': np.array([10.0]),
            'porb': np.array([1000.0]),     # days -> several AU, well detached
            'ecc': np.array([0.1]),
            'metallicity': np.array([0.014]),
        }
        self.assertFalse(default_reject(binary_params)[0])

    def test_low_mass_secondary_rejected(self):
        binary_params = {
            'mass_1': np.array([20.0]),
            'mass_2': np.array([0.02]),     # below the 0.08 Msun minimum
            'porb': np.array([1000.0]),
            'ecc': np.array([0.1]),
            'metallicity': np.array([0.014]),
        }
        self.assertTrue(default_reject(binary_params)[0])

    def test_contact_binary_rejected(self):
        # A 0.5-day orbit puts two massive stars in contact at ZAMS.
        binary_params = {
            'mass_1': np.array([30.0]),
            'mass_2': np.array([25.0]),
            'porb': np.array([0.5]),
            'ecc': np.array([0.0]),
            'metallicity': np.array([0.014]),
        }
        self.assertTrue(default_reject(binary_params)[0])


# ==========================================================================
# Binary-parameter model (sampled + derived coverage)
# ==========================================================================
class TestBinaryModel(unittest.TestCase):

    def _make(self, params, derive_params=None):
        return AdaptiveSampler(
            parameter_space=params, total_systems=10, batch_size=5, BSEDict={},
            is_interesting=lambda bpp: (0, np.array([], dtype=int)),
            derive_params=derive_params, reject_systems=None,
        )

    def test_missing_required_params_raise_at_construction(self):
        # Only porb is sampled; mass_1/mass_2/ecc/metallicity are undefined.
        params = ParameterSpace([Parameter('porb', 10**(0.15), 10**(5.5), dist='sana')])
        with self.assertRaisesRegex(ValueError, 'mass_1'):
            self._make(params)

    def test_derive_params_fills_missing_with_scalars(self):
        # Sample only porb; fix the rest via scalar returns (broadcast to N).
        params = ParameterSpace([Parameter('porb', 10**(0.15), 10**(5.5), dist='sana')])
        sampler = self._make(params, derive_params=lambda s: {
            'mass_1': 30.0, 'mass_2': 25.0, 'ecc': 0.0, 'metallicity': 0.02,
        })
        phys = params.to_physical(params.sample(4, rng=np.random.default_rng(0))[0])
        bp = sampler._binary_params(phys)
        for key in AdaptiveSampler.REQUIRED_PARAMS:
            self.assertEqual(bp[key].shape, (4,))
        np.testing.assert_array_equal(bp['mass_1'], np.full(4, 30.0))

    def test_all_required_sampled_needs_no_derive(self):
        params = ParameterSpace([
            Parameter('mass_1', 5.0, 150.0, dist='kroupa'),
            Parameter('mass_2', 1.0, 100.0, dist='uniform'),
            Parameter('porb', 10**(0.15), 10**(5.5), dist='sana'),
            Parameter('ecc', 1e-9, 0.99, dist='sana_ecc'),
            Parameter('metallicity', 1e-4, 0.03, dist='flat_in_log'),
        ])
        self._make(params)   # must not raise

    def test_derive_params_wrong_length_raises(self):
        params = ParameterSpace([Parameter('porb', 10**(0.15), 10**(5.5), dist='sana')])
        sampler_factory = lambda: self._make(params, derive_params=lambda s: {
            'mass_1': np.array([30.0]),   # wrong length vs the 2-row probe
            'mass_2': 25.0, 'ecc': 0.0, 'metallicity': 0.02,
        })
        with self.assertRaisesRegex(ValueError, 'length'):
            sampler_factory()


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
