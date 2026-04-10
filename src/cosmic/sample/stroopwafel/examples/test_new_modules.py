"""Tests for the new vectorized STROOPWAFEL modules.

Run with: python -m pytest tests/test_new_modules.py -v
Or just: python tests/test_new_modules.py
"""
import sys
import os
import math
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))  # adds COSMIC-stroopwafel/

from stroopwafel import ParameterSpace, Parameter
from stroopwafel.samplers import SAMPLERS
from stroopwafel.priors import PRIORS
from stroopwafel.transforms import to_sampling_space, to_physical_space, transform_bounds
from stroopwafel.mixture_model import GaussianMixture
from stroopwafel.rejection import get_zams_radius, calculate_roche_lobe_radius, default_reject
from stroopwafel.result import STROOPWAFELResult
from stroopwafel.constants import (
    R_COEFF, ZSOL, R_SOL_TO_AU, ALPHA_IMF, SANA_G, SANA_ECC
)


def make_default_params():
    return ParameterSpace([
        Parameter('mass_1', 5.0, 150.0, sampler='kroupa', prior='kroupa'),
        Parameter('q', 0.0, 1.0, sampler='uniform', prior='uniform'),
        Parameter('porb', 0.15, 5.5, sampler='sana', prior='sana'),
        Parameter('ecc', 1e-9, 0.99999999, sampler='sana_ecc', prior='sana_ecc'),
        Parameter('metallicity', 0.0001, 0.03, sampler='flat_in_log', prior='flat_in_log'),
    ])


# ====================================================================
# ParameterSpace tests
# ====================================================================

def test_parameter_space_names_sorted():
    params = make_default_params()
    assert params.names == sorted(params.names)
    assert params.ndim == 5


def test_sampling_produces_valid_arrays():
    params = make_default_params()
    rng = np.random.default_rng(42)
    samples, mask = params.sample(1000, rng=rng)
    assert samples.shape == (1000, 5)
    assert mask.shape == (1000,)
    assert mask.dtype == bool


def test_roundtrip_transform():
    params = make_default_params()
    rng = np.random.default_rng(42)
    samples, _ = params.sample(1000, rng=rng)
    physical = params.to_physical(samples)
    back = params.to_sampling(physical)
    assert np.allclose(samples, back, atol=1e-12)


def test_prior_positive():
    params = make_default_params()
    rng = np.random.default_rng(42)
    samples, mask = params.sample(1000, rng=rng)
    priors = params.compute_prior(samples[mask])
    assert np.all(priors > 0)
    assert np.all(np.isfinite(priors))


def test_in_bounds():
    params = make_default_params()
    rng = np.random.default_rng(42)
    samples, mask = params.sample(1000, rng=rng)
    mask2 = params.in_bounds(samples)
    np.testing.assert_array_equal(mask, mask2)


# ====================================================================
# ZAMS radius / Roche lobe tests (vs old scalar implementation)
# ====================================================================

def old_get_zams_radius(mass, metallicity):
    """Old scalar implementation for comparison."""
    metallicity_xi = math.log10(metallicity / ZSOL)
    rc = []
    for coeff in R_COEFF:
        value = 1; total = 0
        for series in coeff:
            total += series * value
            value *= metallicity_xi
        rc.append(total)
    top = (rc[0] * pow(mass, 2.5) + rc[1] * pow(mass, 6.5)
           + rc[2] * pow(mass, 11) + rc[3] * pow(mass, 19)
           + rc[4] * pow(mass, 19.5))
    bottom = (rc[5] + rc[6] * pow(mass, 2) + rc[7] * pow(mass, 8.5)
              + pow(mass, 18.5) + rc[8] * pow(mass, 19.5))
    return (top / bottom) * R_SOL_TO_AU


def old_roche(m1, m2):
    q = m1 / m2
    return 0.49 / (0.6 + pow(q, -2.0 / 3.0) * math.log(1.0 + pow(q, 1.0 / 3.0)))


def test_zams_radius_matches_old():
    masses = np.linspace(1, 100, 50)
    mets = np.full(50, 0.014)
    new = get_zams_radius(masses, mets)
    old = np.array([old_get_zams_radius(m, z) for m, z in zip(masses, mets)])
    np.testing.assert_allclose(new, old, atol=1e-12)


def test_zams_radius_multiple_metallicities():
    masses = np.array([1.0, 10.0, 50.0])
    mets = np.array([0.001, 0.014, 0.03])
    new = get_zams_radius(masses, mets)
    old = np.array([old_get_zams_radius(m, z) for m, z in zip(masses, mets)])
    np.testing.assert_allclose(new, old, atol=1e-12)


def test_roche_lobe_matches_old():
    m1 = np.array([5.0, 10.0, 20.0, 50.0, 100.0])
    m2 = np.array([3.0, 8.0, 10.0, 25.0, 50.0])
    new = calculate_roche_lobe_radius(m1, m2)
    old = np.array([old_roche(a, b) for a, b in zip(m1, m2)])
    np.testing.assert_allclose(new, old, atol=1e-12)


# ====================================================================
# Default rejection function tests
# ====================================================================

def test_default_reject_basic():
    param_names = ['ecc', 'mass_1', 'metallicity', 'porb', 'q']
    # Create a sample that should NOT be rejected (wide binary)
    samples = np.array([[0.1, 20.0, 0.014, 100.0, 0.5]])  # reasonable binary
    derived = {
        'mass_2': np.array([10.0]),
        'metallicity_1': np.array([0.014]),
        'metallicity_2': np.array([0.014]),
        'separation': np.array([50.0]),  # AU
    }
    rejected = default_reject(samples, derived, param_names)
    assert not rejected[0], "Wide binary should not be rejected"


def test_default_reject_low_mass():
    param_names = ['ecc', 'mass_1', 'metallicity', 'porb', 'q']
    samples = np.array([[0.1, 20.0, 0.014, 100.0, 0.001]])
    derived = {
        'mass_2': np.array([0.02]),  # below minimum secondary mass
        'metallicity_1': np.array([0.014]),
        'metallicity_2': np.array([0.014]),
        'separation': np.array([50.0]),
    }
    rejected = default_reject(samples, derived, param_names)
    assert rejected[0], "Low mass secondary should be rejected"


# ====================================================================
# GaussianMixture tests
# ====================================================================

def test_mixture_from_hits():
    params = make_default_params()
    rng = np.random.default_rng(42)
    samples, _ = params.sample(1000, rng=rng)
    # Use moderate samples as "hits"
    hits = samples[400:410]  # 10 hits

    avg_density = 1.0 / np.power(1000, 1.0 / params.ndim)
    mixture = GaussianMixture.from_hits(hits, params, avg_density)

    assert mixture.n_components == 10
    assert mixture.means.shape == (10, 5)
    assert mixture.covariances.shape == (10, 5, 5)
    assert np.allclose(np.sum(mixture.alphas), 1.0)
    assert np.allclose(mixture.alphas, 0.1)


def test_mixture_pdf_nonnegative():
    params = make_default_params()
    rng = np.random.default_rng(42)
    samples, _ = params.sample(1000, rng=rng)
    hits = samples[400:410]

    avg_density = 1.0 / np.power(1000, 1.0 / params.ndim)
    mixture = GaussianMixture.from_hits(hits, params, avg_density)

    pdf_vals = mixture.pdf(samples[:100])
    assert np.all(pdf_vals >= 0)
    assert np.all(np.isfinite(pdf_vals))


def test_mixture_sample():
    params = make_default_params()
    rng = np.random.default_rng(42)
    samples, _ = params.sample(1000, rng=rng)
    hits = samples[400:410]

    avg_density = 1.0 / np.power(1000, 1.0 / params.ndim)
    mixture = GaussianMixture.from_hits(hits, params, avg_density)

    msamp, mmask, midx = mixture.sample(2000, params, rng=rng)
    assert msamp.shape[1] == 5
    assert len(mmask) == len(msamp)
    assert len(midx) == len(msamp)
    assert np.all(midx >= 0)
    assert np.all(midx < 10)


# ====================================================================
# STROOPWAFELResult tests
# ====================================================================

def test_result_hit_rate():
    result = STROOPWAFELResult()
    result.weights = np.ones(100)
    result.is_hit = np.zeros(100, dtype=bool)
    result.is_hit[:10] = True
    assert abs(result.hit_rate - 0.1) < 1e-10


# ====================================================================
# Prior comparison with old code
# ====================================================================

def test_kroupa_prior_matches():
    """Test that vectorized kroupa prior matches old scalar version."""
    from stroopwafel.priors import kroupa
    lo, hi = 5.0, 150.0
    values = np.linspace(5.1, 149.9, 100)
    vec_result = kroupa(values, lo, hi)

    # Old scalar
    norm = (ALPHA_IMF + 1) / (hi**(ALPHA_IMF + 1) - lo**(ALPHA_IMF + 1))
    old_result = norm * values**ALPHA_IMF
    np.testing.assert_allclose(vec_result, old_result, atol=1e-12)


def test_sana_prior_matches():
    from stroopwafel.priors import sana
    lo, hi = 0.15, 5.5
    values = np.linspace(0.2, 5.4, 100)
    vec_result = sana(values, lo, hi)

    norm = (SANA_G + 1) / (hi**(SANA_G + 1) - lo**(SANA_G + 1))
    old_result = norm * values**SANA_G
    np.testing.assert_allclose(vec_result, old_result, atol=1e-12)


# ====================================================================
# Performance benchmark
# ====================================================================

def bench_zams_radius():
    """Benchmark vectorized vs scalar ZAMS radius."""
    import time
    N = 10000
    masses = np.random.uniform(1, 100, N)
    mets = np.random.uniform(0.0001, 0.03, N)

    # Vectorized
    start = time.time()
    _ = get_zams_radius(masses, mets)
    vec_time = time.time() - start

    # Scalar
    start = time.time()
    _ = [old_get_zams_radius(m, z) for m, z in zip(masses, mets)]
    scalar_time = time.time() - start

    speedup = scalar_time / vec_time if vec_time > 0 else float('inf')
    print(f"\nZAMS radius benchmark (N={N}):")
    print(f"  Vectorized: {vec_time*1000:.1f}ms")
    print(f"  Scalar:     {scalar_time*1000:.1f}ms")
    print(f"  Speedup:    {speedup:.0f}x")


if __name__ == '__main__':
    # Run all tests
    test_funcs = [v for k, v in sorted(globals().items()) if k.startswith('test_')]
    passed = 0
    failed = 0
    for func in test_funcs:
        try:
            func()
            print(f"  [PASS] {func.__name__}")
            passed += 1
        except Exception as e:
            print(f"  [FAIL] {func.__name__}: {e}")
            failed += 1

    print(f"\n{passed} passed, {failed} failed")

    if failed == 0:
        bench_zams_radius()
