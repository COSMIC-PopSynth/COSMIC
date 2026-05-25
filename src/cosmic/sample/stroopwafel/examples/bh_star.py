#!/usr/bin/env python
"""Example: black hole + normal-star binaries with STROOPWAFEL vs Monte Carlo.

A "hit" is any binary where one component is a black hole (kstar=14) while
its companion is still a normal (non-degenerate) star (kstar < 10) and the
system remains bound (sep > 0).  BH+star binaries are intrinsically rare
because they require a massive primary that forms a BH without disrupting
the binary, making them a good testbed for adaptive importance sampling.

The script runs both methods on identical parameter spaces and system budgets,
then prints a side-by-side comparison of hit rates, uncertainties, and
wall-clock times.  The key figure of merit is how many Monte Carlo systems
would be needed to match STROOPWAFEL's statistical precision.

Usage
-----
    python bh_star.py                      # full comparison (default 10 000 systems)
    python bh_star.py --num_systems 50000  # larger run for a clearer signal
    python bh_star.py --mc_only            # MC baseline only
    python bh_star.py --sw_only            # STROOPWAFEL only
"""
import os
import sys
import time
import argparse
import numpy as np

sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
sys.path.append(os.path.join(os.path.dirname(__file__), '../..'))

from stroopwafel import AdaptiveSampler, ParameterSpace, Parameter
from stroopwafel.rejection import default_reject
from stroopwafel import io as swio

# ------------------------------------------------------------------
# CLI
# ------------------------------------------------------------------
parser = argparse.ArgumentParser(
    description=__doc__,
    formatter_class=argparse.RawDescriptionHelpFormatter,
)
parser.add_argument('--num_systems', type=int, default=50000,
                    help='Total systems to evolve per method (default: 50000)')
parser.add_argument('--batch_size', type=int, default=1000,
                    help='Systems per COSMIC call (default: 1000)')
parser.add_argument('--num_cores', type=int, default=8,
                    help='CPU cores for COSMIC (default: 8)')
parser.add_argument('--n_generations', type=int, default=1,
                    help='STROOPWAFEL refinement generations (default: 1)')
parser.add_argument('--mc_only', action='store_true',
                    help='Run Monte Carlo baseline only, skip STROOPWAFEL')
parser.add_argument('--sw_only', action='store_true',
                    help='Run STROOPWAFEL only, skip Monte Carlo baseline')
parser.add_argument('--output_dir', default='output/bh_star',
                    help='Directory for output files (default: output/bh_star)')
parser.add_argument('--seed', type=int, default=117,
                    help='Base random seed; MC uses seed, STROOPWAFEL uses seed+1')
args = parser.parse_args()

# ------------------------------------------------------------------
# BSE physics
# ------------------------------------------------------------------
BSEDict = {
    'xi': 1.0, 'bhflag': 1, 'neta': 0.5, 'windflag': 3, 'wdflag': 1,
    'alpha1': [1.0, 1.0], 'pts1': 0.001, 'pts3': 0.02, 'pts2': 0.01,
    'epsnov': 0.001, 'hewind': 0.5, 'ck': 1000, 'bwind': 0.0,
    'lambdaf': 0.0, 'mxns': 3.0, 'beta': -1.0, 'tflag': 1, 'acc2': 1.5,
    'grflag': 1, 'remnantflag': 4, 'ceflag': 0, 'eddfac': 1.0,
    'ifflag': 0, 'bconst': 3000, 'sigma': 265.0, 'gamma': -2.0,
    'pisn': 45.0,
    # natal_kick_array is intentionally absent: kick parameters are sampled
    # per-binary in the ParameterSpace below, and the engine injects them
    # directly into the InitialBinaryTable for each COSMIC call.
    'bhsigmafrac': 1.0, 'polar_kick_angle': 90,
    'qcrit_array': [0.0] * 16,
    'cekickflag': 2, 'cehestarflag': 0, 'cemergeflag': 0,
    'ecsn': 2.5, 'ecsn_mlow': 1.8, 'aic': 1, 'ussn': 0,
    'sigmadiv': -20.0, 'qcflag': 5, 'eddlimflag': 0,
    'fprimc_array': [2.0 / 21.0] * 16,
    'bhspinflag': 0, 'bhspinmag': 0.0, 'rejuv_fac': 1.0,
    'rejuvflag': 0, 'htpmb': 1, 'ST_cr': 1, 'ST_tide': 1,
    'bdecayfac': 1, 'rembar_massloss': 0.5, 'kickflag': 5,
    'zsun': 0.014, 'bhms_coll_flag': 0, 'don_lim': -1,
    'acc_lim': [-1, -1], 'rtmsflag': 0, 'wd_mass_lim': 1,
    "ppi_co_shift": 0.0, "ppi_extra_ml": 0.0, "fryer_mass_limit": 0,
    "maltsev_mode": 0, "maltsev_fallback": 0.5, "maltsev_pf_prob": 0.1,
    "mm_mu_ns": 800, "mm_mu_bh": 400, "LBV_flag": 1,
    "fryer_fmix": 0.5, "fryer_mcrit_nsbh": 5.0, "smt_periastron_check": 0
}

# ------------------------------------------------------------------
# Parameter space
#
# Orbital / stellar parameters (5 dimensions)
# mass_1       : primary mass [Msun], Kroupa IMF
# q            : mass ratio m2/m1 ∈ [0.01, 1], uniform
# porb         : log10(orbital period / days), Sana power law
#                bounds 0.15 to 5.5 → periods ~1.4 d to ~316 000 d
# ecc          : eccentricity, Sana power law
# metallicity  : metallicity, log-uniform
#
# Primary natal kick magnitude (1 dimension)
# natal_kick_1 : kick speed [km/s], log-normal (mu=5.67, sigma=0.59 in ln-space)
#                physical bounds [0.1, 5000] km/s; median ≈ 291 km/s
#
# Total: 6-dimensional parameter space.
#
# Kick angles (phi_1, theta_1, mean_anomaly_1) and the secondary kick are
# intentionally excluded.  Angles have flat hit-probability across their full
# range so they contribute no information to the mixture model while each
# extra dimension widens the Gaussians by N^(1/D_old - 1/D_new).  The engine
# fills all omitted kick columns with the -100 sentinel so COSMIC draws those
# components from its own prescription (kickflag=5 / sigma=265 km/s).
#
# Note: ParameterSpace sorts parameters alphabetically, so the internal
# column order is fixed and independent of the order given here.
# ------------------------------------------------------------------
params = ParameterSpace([
    # --- orbital / stellar ---
    Parameter('mass_1',       5.0,    150.0,      sampler='kroupa',      prior='kroupa'),
    Parameter('q',            0.01,   1.0,        sampler='uniform',     prior='uniform'),
    Parameter('porb',         0.15,   5.5,        sampler='sana',        prior='sana'),
    Parameter('ecc',          1e-9,   0.99999999, sampler='sana_ecc',    prior='sana_ecc'),
    Parameter('metallicity',  0.0001, 0.03,       sampler='flat_in_log', prior='flat_in_log'),
    # --- primary natal kick magnitude only ---
    # Parameter('natal_kick_1', 0.1,    100.0,     sampler='log_normal',  prior='log_normal'),
])

# ------------------------------------------------------------------
# Derived quantities
# ------------------------------------------------------------------
def compute_derived(samples_physical, param_names):
    """Compute mass_2, metallicities, and orbital separation from samples.

    Parameters
    ----------
    samples_physical : numpy.ndarray
        (N, D) array of samples in physical space.
    param_names : list of str
        Column labels for each dimension.

    Returns
    -------
    dict
        Keys: 'mass_2', 'metallicity_1', 'metallicity_2', 'separation'.
    """
    idx = {name: i for i, name in enumerate(param_names)}
    m1   = samples_physical[:, idx['mass_1']]
    q    = samples_physical[:, idx['q']]
    porb = samples_physical[:, idx['porb']]   # days (sana sampler output)
    z    = samples_physical[:, idx['metallicity']]

    mass_2 = m1 * q

    # Kepler's third law: a³ [AU³] = (P [yr])² · M [Msun]
    # Convert period from days to years before applying.
    separation = ((porb / 365.25) ** 2 * (m1 + mass_2)) ** (1.0 / 3.0)

    return {
        'mass_2':        mass_2,
        'metallicity_1': z,
        'metallicity_2': z,
        'separation':    separation,
    }

# ------------------------------------------------------------------
# Hit definition: BH (kstar=14) + normal star (kstar 0–9), still bound
# ------------------------------------------------------------------
_STAR_KSTARS = set(range(10))   # kstar 0–9: non-degenerate stars

def is_bh_star(bpp):
    """Identify binaries that are in a bound BH + normal-star phase for ≥ 100 Myr.

    Parameters
    ----------
    bpp : pandas.DataFrame
        COSMIC binary population parameters output.

    Returns
    -------
    n_hits : int
        Number of distinct binaries satisfying the criterion.
    hit_bin_nums : numpy.ndarray
        Integer bin_num values of those binaries.
    """
    bh_star_mask = (
        (   (bpp.kstar_1 == 14) & bpp.kstar_2.isin(_STAR_KSTARS))
        | (bpp.kstar_1.isin(_STAR_KSTARS) & (bpp.kstar_2 == 14))
    ) & (bpp.sep > 0)

    bh_star_rows = bpp.loc[bh_star_mask]
    if len(bh_star_rows) == 0:
        return 0, np.array([], dtype=int)

    # Group by bin_num (the natural key) so that min/max tphys are aligned
    # on the same index.  drop_duplicates would give rows with different
    # integer-row indices that pandas would NOT align correctly on subtraction.
    phase_start = bh_star_rows.groupby('bin_num')['tphys'].min()  # Series indexed by bin_num
    phase_end   = bh_star_rows.groupby('bin_num')['tphys'].max()  # Series indexed by bin_num
    duration    = phase_end - phase_start                          # tphys in Myr

    long_enough_bin_nums = duration[duration >= 100.0].index.values
    return len(long_enough_bin_nums), long_enough_bin_nums

# ------------------------------------------------------------------
# Helper: run one sampler and return (result, elapsed_seconds)
# ------------------------------------------------------------------
def run_sampler(mc_only, seed):
    label = "Monte Carlo" if mc_only else "STROOPWAFEL"
    print(f"\n{'='*60}")
    print(f"  Running {label}  (seed={seed})")
    print(f"{'='*60}")

    sw = AdaptiveSampler(
        parameter_space=params,
        total_systems=args.num_systems,
        batch_size=args.batch_size,
        BSEDict=BSEDict,
        compute_derived=compute_derived,
        reject_systems=default_reject,
        is_interesting=is_bh_star,
        output_path=os.path.join(args.output_dir, 'mc' if mc_only else 'sw'),
        nproc=args.num_cores,
        n_generations=args.n_generations,
        mc_only=mc_only,
        seed=seed,
    )

    t0 = time.time()
    result = sw.run()
    elapsed = time.time() - t0

    return result, elapsed

# ------------------------------------------------------------------
# Helper: summarise one result
# ------------------------------------------------------------------
def print_summary(label, result, elapsed):
    raw_hits = int(np.sum(result.is_hit))
    hit_rate = result.hit_rate
    uncertainty = result.hit_rate_uncertainty

    print(f"\n--- {label} summary ---")
    print(f"  Systems evolved  : {len(result.weights):,}")
    print(f"  Raw hits found   : {raw_hits:,}")
    print(f"  Weighted hit rate: {hit_rate:.6e} ± {uncertainty:.6e}")
    print(f"  Wall-clock time  : {elapsed:.1f} s")

# ------------------------------------------------------------------
# Helper: side-by-side comparison
# ------------------------------------------------------------------
def print_comparison(mc_result, mc_elapsed, sw_result, sw_elapsed):
    mc_rate  = mc_result.hit_rate
    mc_unc   = mc_result.hit_rate_uncertainty
    sw_rate  = sw_result.hit_rate
    sw_unc   = sw_result.hit_rate_uncertainty

    print(f"\n{'='*60}")
    print("  Efficiency comparison")
    print(f"{'='*60}")
    print(f"{'':30s}  {'Monte Carlo':>15s}  {'STROOPWAFEL':>15s}")
    print(f"  {'Systems evolved':<28s}  {len(mc_result.weights):>15,}  {len(sw_result.weights):>15,}")
    print(f"  {'Raw hits found':<28s}  {int(np.sum(mc_result.is_hit)):>15,}  {int(np.sum(sw_result.is_hit)):>15,}")
    print(f"  {'Weighted hit rate':<28s}  {mc_rate:>15.4e}  {sw_rate:>15.4e}")
    print(f"  {'Uncertainty (1σ)':<28s}  {mc_unc:>15.4e}  {sw_unc:>15.4e}")
    print(f"  {'Wall-clock time (s)':<28s}  {mc_elapsed:>15.1f}  {sw_elapsed:>15.1f}")

    if sw_unc > 0 and mc_unc > 0:
        # How many MC systems would give the same precision as SW?
        equivalent_mc = len(mc_result.weights) * (mc_unc / sw_unc) ** 2
        speedup = equivalent_mc / len(sw_result.weights)
        print(f"\n  STROOPWAFEL is ~{speedup:.1f}x more statistically efficient:")
        print(f"  Monte Carlo would need ~{equivalent_mc:,.0f} systems to match")
        print(f"  STROOPWAFEL's precision on {len(sw_result.weights):,} systems.")

# ------------------------------------------------------------------
# Main
# ------------------------------------------------------------------
if __name__ == '__main__':
    os.makedirs(args.output_dir, exist_ok=True)

    mc_result = mc_elapsed = None
    sw_result = sw_elapsed = None

    if not args.sw_only:
        mc_result, mc_elapsed = run_sampler(mc_only=True,  seed=args.seed)
        print_summary("Monte Carlo", mc_result, mc_elapsed)
        swio.save_result(os.path.join(args.output_dir, 'mc_result.h5'), mc_result)

    if not args.mc_only:
        sw_result, sw_elapsed = run_sampler(mc_only=False, seed=args.seed + 1)
        print_summary("STROOPWAFEL", sw_result, sw_elapsed)
        swio.save_result(os.path.join(args.output_dir, 'sw_result.h5'), sw_result)

    if mc_result is not None and sw_result is not None:
        print_comparison(mc_result, mc_elapsed, sw_result, sw_elapsed)
