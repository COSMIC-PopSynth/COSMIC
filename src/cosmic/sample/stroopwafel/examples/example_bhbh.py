#!/usr/bin/env python
"""Example: Find merging BH-BH binaries using the new vectorized STROOPWAFEL API.

Equivalent to the old tests/BHBH_testing.py but using the new API.
"""
import os
import sys
import time
import argparse
import numpy as np

sys.path.append(os.path.join(os.path.dirname(__file__), '..'))  # adds COSMIC-stroopwafel/

from stroopwafel import AdaptiveSampler, ParameterSpace, Parameter
from stroopwafel.presets import merging_dco
from stroopwafel.rejection import default_reject
from stroopwafel import io as swio

# ------------------------------------------------------------------
# Parse arguments
# ------------------------------------------------------------------
parser = argparse.ArgumentParser()
parser.add_argument('--num_systems', type=int, default=10000)
parser.add_argument('--num_cores', type=int, default=1)
parser.add_argument('--num_per_core', type=int, default=100)
parser.add_argument('--mc_only', type=bool, default=False)
parser.add_argument('--output_dir', default='output/BHBH_new')
parser.add_argument('--model', default='fiducial')
parser.add_argument('--seed', type=int, default=None)
args = parser.parse_args()

# ------------------------------------------------------------------
# Physics models (same as old code)
# ------------------------------------------------------------------
fiducial = {
    'xi': 1.0, 'bhflag': 1, 'neta': 0.5, 'windflag': 3, 'wdflag': 1,
    'alpha1': 1.0, 'pts1': 0.001, 'pts3': 0.02, 'pts2': 0.01,
    'epsnov': 0.001, 'hewind': 0.5, 'ck': 1000, 'bwind': 0.0,
    'lambdaf': 0.0, 'mxns': 3.0, 'beta': -1.0, 'tflag': 1, 'acc2': 1.5,
    'grflag': 1, 'remnantflag': 4, 'ceflag': 0, 'eddfac': 1.0,
    'ifflag': 0, 'bconst': 3000, 'sigma': 265.0, 'gamma': -2.0,
    'pisn': 45.0,
    'natal_kick_array': [[-100.0, -100.0, -100.0, -100.0, 0.0],
                         [-100.0, -100.0, -100.0, -100.0, 0.0]],
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
    'acc_lim': -1, 'rtmsflag': 0, 'wd_mass_lim': 1,
}

BSEDict = fiducial  # extend with model variants as needed

# ------------------------------------------------------------------
# Define parameter space
# ------------------------------------------------------------------
params = ParameterSpace([
    Parameter('mass_1', 5.0, 150.0, sampler='kroupa', prior='kroupa'),
    Parameter('q', 0.0, 1.0, sampler='uniform', prior='uniform'),
    Parameter('porb', 0.15, 5.5, sampler='sana', prior='sana'),
    Parameter('ecc', 1e-9, 0.99999999, sampler='sana_ecc', prior='sana_ecc'),
    Parameter('metallicity', 0.0001, 0.03, sampler='flat_in_log', prior='flat_in_log'),
])

# ------------------------------------------------------------------
# Define derived quantities (vectorized)
# ------------------------------------------------------------------
def compute_derived(samples_physical, param_names):
    """Compute derived quantities from sampled parameters."""
    idx = {name: i for i, name in enumerate(param_names)}
    m1 = samples_physical[:, idx['mass_1']]
    q = samples_physical[:, idx['q']]
    porb = samples_physical[:, idx['porb']]
    z = samples_physical[:, idx['metallicity']]

    mass_2 = m1 * q
    separation = ((porb ** 2) * (m1 + mass_2)) ** (1.0 / 3.0)

    return {
        'mass_2': mass_2,
        'metallicity_1': z,
        'metallicity_2': z,
        'separation': separation,
    }

# ------------------------------------------------------------------
# Hit definition: merging BH-BH binaries within Hubble time
# ------------------------------------------------------------------
is_interesting = merging_dco(kstar_1=[14], kstar_2=[14], max_merge_time=13.7)

# ------------------------------------------------------------------
# Run
# ------------------------------------------------------------------
if __name__ == '__main__':
    start = time.time()

    sw = AdaptiveSampler(
        parameter_space=params,
        total_systems=args.num_systems,
        batch_size=args.num_per_core,
        BSEDict=BSEDict,
        compute_derived=compute_derived,
        reject_systems=default_reject,
        is_interesting=is_interesting,
        output_path=args.output_dir,
        nproc=args.num_cores,
        mc_only=args.mc_only,
        seed=args.seed,
    )

    result = sw.run()

    # Save results
    output_file = os.path.join(args.output_dir, 'stroopwafel_result.h5')
    swio.save_result(output_file, result)

    elapsed = time.time() - start
    print(f"\nTotal time: {elapsed:.1f}s")
    print(f"Results saved to {output_file}")
