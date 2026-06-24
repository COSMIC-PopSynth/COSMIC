"""
``windflag``
============

This example shows the effect of the ``windflag`` setting on massive-star
evolution. ``windflag`` selects the stellar wind mass-loss prescription used by
COSMIC. The initial metallicity and binary setup are held fixed, while only
``windflag`` is changed.

The systems are very wide binaries with low-mass companions, so the example is
focused on stellar winds rather than Roche-lobe overflow, common-envelope
evolution, or tides.
"""

import sys

sys.path.append("..")
import generate_default_bsedict

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from cosmic.evolve import Evolve
from cosmic.sample import InitialBinaryTable

BSEDict = generate_default_bsedict.get_default_BSE_settings(to_python=True)
plt.style.use("../_static/gallery.mplstyle")


windflag_labels = {
    -1: "-1 no winds",
    0: "0 SSE/BSE",
    1: "1 StarTrack",
    2: "2 Vink OB/WR",
    3: "3 Vink+LBV",
    5: "5 0.33x flag 3",
    6: "6 Bjorklund+23",
    7: "7 Krticka+24",
}

windflag_values = list(windflag_labels)
masses = np.arange(5.0, 81.0, 1.0)

kstar_labels = {
    1: "1: MS",
    2: "2: HG",
    3: "3: GB",
    4: "4: CHeB",
    5: "5: EAGB",
    6: "6: TPAGB",
    7: "7: HeMS",
    8: "8: HeHG",
    9: "9: HeGB",
    10: "10: HeWD",
    11: "11: COWD",
    12: "12: ONeWD",
    13: "13: NS",
    14: "14: BH",
    15: "15: MR",
}


def make_wide_binary_grid():
    """Build the wide-binary grid used to isolate stellar wind effects."""
    return InitialBinaryTable.InitialBinaries(
        m1=masses,
        m2=np.ones_like(masses) * 0.1,
        porb=np.ones_like(masses) * 1e6,
        ecc=np.zeros_like(masses),
        tphysf=np.ones_like(masses) * 13700.0,
        kstar1=np.ones_like(masses),
        kstar2=np.ones_like(masses),
        metallicity=np.ones_like(masses) * 0.014,
    )


def summarize_final_state(bpp, bcm, windflag):
    """Collect final stellar quantities for one wind prescription."""
    final_bpp = bpp.sort_values(["bin_num", "tphys"]).groupby("bin_num").last()
    final_bcm = bcm.sort_values(["bin_num", "tphys"]).groupby("bin_num").last()

    summary = pd.DataFrame(index=final_bcm.index)
    summary["windflag"] = windflag
    summary["windflag_label"] = windflag_labels[windflag]
    summary["initial_mass_1"] = masses[summary.index.astype(int)]
    summary["final_mass_1"] = final_bcm["mass_1"]
    summary["final_co_core_1"] = final_bcm["massc_co_layer_1"]
    summary["final_kstar_1"] = final_bcm["kstar_1"]
    summary["final_evol_type"] = final_bpp["evol_type"]
    return summary.reset_index(drop=True)


def evolve_windflag_grid(binary_grid):
    """Run the same initial grid for each windflag value."""
    summaries = []

    for windflag in windflag_values:
        settings = BSEDict.copy()
        settings["windflag"] = windflag
        settings["random_seed"] = 1

        bpp, bcm, initC, kick_info = Evolve.evolve(
            initialbinarytable=binary_grid,
            BSEDict=settings,
        )
        summaries.append(summarize_final_state(bpp, bcm, windflag))

    return pd.concat(summaries, ignore_index=True)


binary_grid = make_wide_binary_grid()
summaries = evolve_windflag_grid(binary_grid)

fig, axes = plt.subplots(1, 2, figsize=(13, 5), sharex=True, constrained_layout=True)

for offset, windflag in zip(np.linspace(-0.14, 0.14, len(windflag_values)), windflag_values):
    subset = summaries[summaries["windflag"].eq(windflag)]
    label = windflag_labels[windflag]
    axes[0].plot(subset["initial_mass_1"], subset["final_mass_1"], label=label)
    axes[1].scatter(
        subset["initial_mass_1"],
        subset["final_kstar_1"] + offset,
        s=16,
        label=label,
    )

axes[0].set_xlabel("Initial mass [$M_\\odot$]")
axes[0].set_ylabel("Final mass [$M_\\odot$]")
axes[0].set_title("Final mass")

axes[1].set_xlabel("Initial mass [$M_\\odot$]")
axes[1].set_ylabel("Final kstar")
axes[1].set_title("Final kstar type")

used_kstars = sorted(summaries["final_kstar_1"].astype(int).unique())
axes[1].set_yticks(used_kstars)
axes[1].set_yticklabels([kstar_labels.get(k, str(k)) for k in used_kstars])

for ax in axes:
    ax.legend(title="windflag", fontsize=7, title_fontsize=8, markerscale=1.2)
    ax.grid(True, alpha=0.3)

fig.suptitle("Effect of windflag at fixed Z=0.014 and zsun=0.014")
plt.savefig("windflag_flagtest.png", dpi=150)
plt.show()
