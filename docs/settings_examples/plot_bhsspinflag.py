"""
bhsspinflag/bhspinmag
=======================
These plots test ``bhspinflag`` and ``bhspinmag``, which set the model used to assign a spin to every 
black hole remnant. 

``bhspinflag=0`` gives every BH a fixed spin equal to ``bhspinmag``. 
``bhspinflag=1`` gives each BH's spin uniformly at random between 0 and ``bhspinmag``. 
``bhsspinflag=2`` sets spin from the progenitor's corse mass instead, so it decreases with initial mass and does use ``bhspinmag``. 

A grid of single stars (ZAMS mass 20-100 Msun, windflag=-1, remnantflag=4,
Z=0.002) is evolved once per ``bhspinmag`` value in {0, 0.2, 0.5, 0.8, 1.0} via
``get_bh_spins``, which reads each star's final ``bhspin_1`` once it has become
a BH (``kstar_1 == 14``). The left panel shows the fixed-spin values under
``bhspinflag=0``; the middle panel shows histograms of the random spins under
``bhspinflag=1``; the right panel shows the mass-dependent spin curve under
``bhspinflag=2``, which is independent of ``bhspinmag``.
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from cosmic.sample.initialbinarytable import InitialBinaryTable
from cosmic.evolve import Evolve

n_grid = 100  
masses = np.linspace(20.0, 100.0, n_grid)

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
    "mxns": 3.0, "fryer_fmix": 1.0,
    "fryer_mcrit_nsbh": 5.75, "rembar_massloss": 0.5, "wd_mass_lim": 1,
    "maltsev_mode": 0, "maltsev_fallback": 0.5, "maltsev_pf_prob": 0.1,
    "pisn": -2, "ppi_co_shift": 0.0, "ppi_extra_ml": 0.0,
    "bhspinflag": 0, "bhspinmag": 0.0,
    "grflag": 1, "eddfac": 10, "gamma": -2, "don_lim": -1,
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

def get_bh_spins(BSEDict, n_grid):
    bpp, bcm, initC, kick_info = Evolve.evolve(
        initialbinarytable=binary_grid, BSEDict=BSEDict, SSEDict=SSEDict
    )
    init_masses, spins = [], []
    for i in range(n_grid):
        rows = bpp.loc[i]
        init_mass = rows.iloc[0]['mass_1']
        bh_rows = rows[rows['kstar_1'] == 14]
        if len(bh_rows) > 0:
            spin = bh_rows.iloc[-1]['bhspin_1']
            init_masses.append(init_mass)
            spins.append(spin)
    return init_masses, spins

spinmag_values = [0.0, 0.2, 0.5, 0.8, 1.0]
colors = ['red', 'orange', 'green', 'blue', 'purple']

fig, axes = plt.subplots(1, 3, figsize=(18, 6))

# ── panel 1: bhspinflag=0 — bar chart showing fixed spin values ──────────────
fixed_spins = []
for spinmag in spinmag_values:
    BSEDict = BSEDict_base.copy()
    BSEDict['bhspinflag'] = 0
    BSEDict['bhspinmag'] = spinmag
    init_masses, spins = get_bh_spins(BSEDict, n_grid)
    if len(spins) > 0:
        fixed_spins.append(spins[0])  # all spins are the same so just take first
    print(f"bhspinflag=0, bhspinmag={spinmag}: spin={spins[0]:.2f}")

bars = axes[0].bar(
    [str(s) for s in spinmag_values],
    fixed_spins,
    color=colors,
    edgecolor='black'
)
# add value labels on top of bars
for bar, spin in zip(bars, fixed_spins):
    axes[0].text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.01,
                f'{spin:.1f}', ha='center', va='bottom', fontsize=11)

axes[0].set_xlabel('bhspinmag', fontsize=11)
axes[0].set_ylabel('BH Spin', fontsize=11)
axes[0].set_title('bhspinflag=0', fontsize=13)
axes[0].set_ylim(0, 1.15)
axes[0].grid(True, alpha=0.3, axis='y')

# ── panel 2: bhspinflag=1 — histogram of random spins ────────────────────────
for spinmag, color in zip(spinmag_values[1:], colors[1:]):  # skip 0.0
    BSEDict = BSEDict_base.copy()
    BSEDict['bhspinflag'] = 1
    BSEDict['bhspinmag'] = spinmag
    init_masses, spins = get_bh_spins(BSEDict, n_grid)
    if len(spins) > 0:
        axes[1].hist(spins, bins=10, alpha=1.0, color=color,
                    label=f'bhspinmag={spinmag}', range=(0, 1))
        print(f"bhspinflag=1, bhspinmag={spinmag}: "
              f"spin range={min(spins):.2f}-{max(spins):.2f}")

axes[1].set_xlabel('BH Spin', fontsize=11)
axes[1].set_ylabel('Count', fontsize=11)
axes[1].set_title('bhspinflag=1', fontsize=13)
axes[1].legend(fontsize=11)
axes[1].grid(True, alpha=0.3)

# ── panel 3: bhspinflag=2 — core mass dependent spin ─────────────────────────
BSEDict = BSEDict_base.copy()
BSEDict['bhspinflag'] = 2
BSEDict['bhspinmag'] = 1.0
init_masses, spins = get_bh_spins(BSEDict, n_grid)
if len(init_masses) > 0:
    axes[2].plot(init_masses, spins, 'o-', color='steelblue',
                label='Core mass dependent spin')
    print(f"bhspinflag=2: {len(spins)} BHs, "
          f"spin range={min(spins):.2f}-{max(spins):.2f}")

axes[2].set_xlabel('Initial Stellar Mass (Msun)', fontsize=11)
axes[2].set_ylabel('BH Spin', fontsize=11)
axes[2].set_title('bhspinflag=2', fontsize=13)
axes[2].set_ylim(-0.05, 1.05)
axes[2].legend(fontsize=11)
axes[2].grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('bhspinflag.png', dpi=150)
plt.show()
