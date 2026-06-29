#!/usr/bin/env python
import time
import argparse

from cosmic.sample.stroopwafel import AdaptiveSampler, ParameterSpace, Parameter
from cosmic.sample.stroopwafel.presets import merging_dco

import sys
sys.path.append("../../../../../docs")
import generate_default_bsedict
BSEDict = generate_default_bsedict.get_default_BSE_settings(to_python=True)


parser = argparse.ArgumentParser()
parser.add_argument('--num_systems', type=int, default=100_000)
parser.add_argument('--num_cores', type=int, default=1)
parser.add_argument('--num_per_core', type=int, default=1000)
parser.add_argument('--mc_only', type=bool, default=False)
parser.add_argument('--seed', type=int, default=None)
args = parser.parse_args()

# COSMIC v4+ single stellar evolution settings (sse engine; swap for METISSE)
SSEDict = {'stellar_engine': 'sse'}

# ------------------------------------------------------------------
# Define parameter space
# ------------------------------------------------------------------
params = ParameterSpace([
    Parameter('mass_1', 5.0, 150.0, dist='kroupa'),
    Parameter('q', 0.0, 1.0, dist='uniform'),
    Parameter('porb', 10**(0.15), 10**(5.5), dist='sana'),
    Parameter('ecc', 1e-9, 0.99999999, dist='sana_ecc'),
    Parameter('metallicity', 0.0001, 0.03, dist='flat_in_log'),
    Parameter('natal_kick_1', 0.0, 5000.0, dist='disberg')
])

def derive_params(sampled):
    """Return the one required parameter (mass_2) not sampled directly."""
    return {'mass_2': sampled['mass_1'] * sampled['q']}

is_interesting = merging_dco(kstar_1=[14], kstar_2=[14], max_merge_time=13.7)

if __name__ == '__main__':
    start = time.time()

    sw = AdaptiveSampler(
        parameter_space=params,
        total_systems=args.num_systems,
        batch_size=args.num_per_core,
        BSEDict=BSEDict,
        SSEDict=SSEDict,
        is_interesting=is_interesting,
        derive_params=derive_params,
        reject_systems="default",
        nproc=args.num_cores,
        mc_only=args.mc_only,
        seed=args.seed,
    )

    result = sw.run()
    result.save("output/bh_bh/result.h5")

    elapsed = time.time() - start
    print(f"\nTotal time: {elapsed:.1f}s")
    print(f"Results saved to output/bh_bh/result.h5")
