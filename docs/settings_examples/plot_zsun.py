"""
``zsun``
========

This example shows the effect of the ``zsun`` setting at fixed absolute
metallicity. ``zsun`` sets the reference value COSMIC uses for solar
metallicity. It does not change the stellar metallicity input directly;
instead, it changes the ratio ``Z / zsun`` used by metallicity-dependent wind
prescriptions.

The stars all have ``Z = 0.014``. Changing ``zsun`` therefore changes how
metal-rich the same stars appear relative to solar, which can change wind mass
loss and the final stellar or remnant type.
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


zsun_values = [0.010, 0.014, 0.019, 0.020, 0.030]
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
    """Build a very wide binary grid to isolate wind/remnant effects."""
    n_systems = len(masses)
    return InitialBinaryTable.InitialBinaries(
        m1=masses,
        m2=np.ones(n_systems) * 0.1,
        porb=np.ones(n_systems) * 1e6,
        ecc=np.zeros(n_systems),
        tphysf=np.ones(n_systems) * 13700.0,
        kstar1=np.ones(n_systems),
        kstar2=np.ones(n_systems),
        metallicity=np.ones(n_systems) * 0.014,
    )


def summarize_final_state(bpp, bcm, zsun):
    """Collect final primary quantities for one zsun value."""
    final_bpp = bpp.sort_values(["bin_num", "tphys"]).groupby("bin_num").last()
    final_bcm = bcm.sort_values(["bin_num", "tphys"]).groupby("bin_num").last()

    summary = pd.DataFrame(index=final_bcm.index)
    summary["zsun"] = zsun
    summary["initial_mass_1"] = masses[summary.index.astype(int)]
    summary["final_mass_1"] = final_bcm["mass_1"]
    summary["mass_lost_1"] = summary["initial_mass_1"] - summary["final_mass_1"]
    summary["final_kstar_1"] = final_bcm["kstar_1"]
    summary["final_evol_type"] = final_bpp["evol_type"]
    return summary.reset_index(drop=True)


def evolve_zsun_grid(binary_grid):
    """Run the same initial grid for each zsun value."""
    summaries = []

    for zsun in zsun_values:
        settings = BSEDict.copy()
        settings["zsun"] = zsun
        settings["random_seed"] = 1

        bpp, bcm, initC, kick_info = Evolve.evolve(
            initialbinarytable=binary_grid,
            BSEDict=settings,
        )
        summaries.append(summarize_final_state(bpp, bcm, zsun))

    return pd.concat(summaries, ignore_index=True)


binary_grid = make_wide_binary_grid()
summaries = evolve_zsun_grid(binary_grid)

fig, axes = plt.subplots(1, 2, figsize=(13, 5), sharex=True, constrained_layout=True)

for offset, zsun in zip(np.linspace(-0.10, 0.10, len(zsun_values)), zsun_values):
    subset = summaries[summaries["zsun"].eq(zsun)].sort_values("initial_mass_1")
    label = f"zsun={zsun:g}"

    axes[0].plot(
        subset["initial_mass_1"],
        subset["final_mass_1"],
        label=label,
    )
    axes[1].scatter(
        subset["initial_mass_1"],
        subset["final_kstar_1"] + offset,
        s=14,
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
    ax.legend(title="zsun", fontsize=8, title_fontsize=9, markerscale=1.2)
    ax.grid(True, alpha=0.3)

fig.suptitle("Effect of zsun at fixed absolute metallicity Z=0.014")
plt.savefig("zsun_flagtest.png", dpi=150)
plt.show()
