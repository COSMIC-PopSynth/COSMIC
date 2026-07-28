"""
mxns
====

This example tests ``mxns``, which sets the highest mass a neutron star is
allowed to have. Any remnant that would end up heavier than ``mxns`` is
turned into a black hole instead of a neutron star.

Raising ``mxns`` lets more of the heavier remnants stay classified as
neutron stars rather than black holes, so a higher ``mxns`` means more
neutron stars and fewer black holes come out of the same underlying
population of stars.

A grid of single stars is evolved once for each of several ``mxns`` values,
from 1.5 up to 7.0 solar masses, and the number of stars ending up as a
neutron star versus a black hole is counted and plotted as a bar chart for
each value.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from cosmic.sample.initialbinarytable import InitialBinaryTable
from cosmic.evolve import Evolve


n_grid = 30
masses = np.linspace(5.0, 50.0, n_grid)

SSEDict = {'stellar_engine': 'sse'}

BSEDict_base = {
    "pts1": 0.001, "pts2": 0.01, "pts3": 0.02, "zsun": 0.014,
    "windflag": -1,
    "eddlimflag": 0, "neta": 0.5, "bwind": 0.0, "hewind": 0.5, "beta": 0.125,
    "xi": 0.5, "acc2": 1.5, "LBV_flag": 1, "alpha1": [1.0, 1.0],
    "lambdaf": 0.0, "ceflag": 1, "cekickflag": 2, "cemergeflag": 1,
    "cehestarflag": 0, "qcflag": 5,
    "qcrit_array": [0.0]*16,
    "kickflag": 5, "sigma": 265.0, "bhflag": 1, "bhsigmafrac": 1.0,
    "sigmadiv": -20.0, "ecsn": 2.25, "ecsn_mlow": 1.6, "aic": 1, "ussn": 1,
    "polar_kick_angle": 90.0,
    "natal_kick_array": [[-100.0,-100.0,-100.0,-100.0,0.0],[-100.0,-100.0,-100.0,-100.0,0.0]],
    "mm_mu_ns": 400.0, "mm_mu_bh": 200.0,
    "remnantflag": 4,
    "fryer_mass_limit": 0,
    "mxns": 3.0,
    "fryer_fmix": 1.0,
    "fryer_mcrit_nsbh": 5.75, "rembar_massloss": 0.5, "wd_mass_lim": 1,
    "maltsev_mode": 0, "maltsev_fallback": 0.5, "maltsev_pf_prob": 0.1,
    "pisn": -2, "ppi_co_shift": 0.0, "ppi_extra_ml": 0.0, "bhspinflag": 0,
    "bhspinmag": 0.0, "grflag": 1, "eddfac": 10, "gamma": -2, "don_lim": -1,
    "acc_lim": [-1, -1], "smt_periastron_check": 0, "tflag": 1, "ST_tide": 1,
    "fprimc_array": [2.0/21.0]*16,
    "ifflag": 1, "wdflag": 1, "epsnov": 0.001, "bdecayfac": 1,
    "bconst": 3000, "ck": 1000, "rejuv_fac": 1.0, "rejuvflag": 0,
    "bhms_coll_flag": 0, "htpmb": 1, "ST_cr": 1, "rtmsflag": 0
}

binary_grid = InitialBinaryTable.InitialBinaries(
    m1=masses,
    m2=np.ones(n_grid) * 0.1,
    porb=np.ones(n_grid) * 100000.0,
    ecc=np.zeros(n_grid),
    tphysf=np.ones(n_grid) * 13700.0,
    kstar1=np.ones(n_grid),
    kstar2=np.ones(n_grid),
    metallicity=np.ones(n_grid) * 0.002
)

mxns_values = [1.5, 2.0, 3.0, 5.0, 7.0]

ns_counts = []
bh_counts = []

for mxns in mxns_values:
    BSEDict = BSEDict_base.copy()
    BSEDict['mxns'] = mxns

    bpp, bcm, initC, kick_info = Evolve.evolve(
        initialbinarytable=binary_grid, BSEDict=BSEDict, SSEDict=SSEDict
    )

    ns_count = 0
    bh_count = 0
    for i in range(n_grid):
        rows = bpp.loc[i]
        final_kstar = rows.iloc[-1]['kstar_1']
        if final_kstar == 13:
            ns_count += 1
        elif final_kstar == 14:
            bh_count += 1

    ns_counts.append(ns_count)
    bh_counts.append(bh_count)
    print(f"mxns={mxns}: {ns_count} NS, {bh_count} BH")

# bar chart
x = np.arange(len(mxns_values))
width = 0.35

fig, ax = plt.subplots(figsize=(10, 6))

bars_ns = ax.bar(x - width/2, ns_counts, width, 
                  label='Neutron Stars', color='steelblue')
bars_bh = ax.bar(x + width/2, bh_counts, width,
                  label='Black Holes', color='firebrick')

for bar in bars_ns:
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.2,
            str(int(bar.get_height())), ha='center', va='bottom', fontsize=11)
for bar in bars_bh:
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.2,
            str(int(bar.get_height())), ha='center', va='bottom', fontsize=11)

ax.set_xlabel('mxns (Msun)', fontsize=12)
ax.set_ylabel('Number of Remnants', fontsize=12)
ax.set_title('Effect of mxns on NS vs BH Count\n(higher mxns = more NS, fewer BH)', fontsize=13)
ax.set_xticks(x)
ax.set_xticklabels([f'mxns={m}\n(default)' if m == 3.0 else f'mxns={m}'
                    for m in mxns_values])
ax.legend(fontsize=11)
ax.grid(True, alpha=0.3, axis='y')

plt.tight_layout()
plt.savefig('mxns.png', dpi=150)
plt.show()
