"""
``hewind``
==========

This example shows the effect of the ``hewind`` setting on naked helium-star
evolution. ``hewind`` scales the wind mass loss for naked helium stars, so this
test initializes the primary as a helium main-sequence star with ``kstar_1 = 7``.

The binaries are very wide, which keeps the example focused on helium-star wind
mass loss rather than Roche-lobe overflow or common-envelope evolution. The
notebook version of this test found the clearest direct ``hewind`` dependence
with ``windflag = 0``, so this script uses that setup for the final diagnostic.
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
plt.rcParams.update({
    "font.size": 9,
    "axes.titlesize": 11,
    "axes.labelsize": 10,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "legend.fontsize": 8,
    "legend.title_fontsize": 8,
})


hewind_values = [0.0, 0.25, 0.5, 0.75, 1.0]
he_masses = np.array([3.0, 5.0, 8.0, 12.0, 20.0])


def hewind_label(hewind):
    """Format labels so the no-hewind and default cases are explicit."""
    if hewind == 0.0:
        return "0 (none)"
    if hewind == 0.5:
        return "0.5 (default)"
    return f"{hewind:g}"


def make_helium_star_grid():
    """Build wide binaries with naked helium-star primaries."""
    n_systems = len(he_masses)
    return InitialBinaryTable.InitialBinaries(
        m1=he_masses,
        m2=np.ones(n_systems),
        porb=np.ones(n_systems) * 1.0e6,
        ecc=np.zeros(n_systems),
        tphysf=np.ones(n_systems) * 2000.0,
        kstar1=np.ones(n_systems) * 7,
        kstar2=np.ones(n_systems),
        metallicity=np.ones(n_systems) * 0.014,
    )


def summarize_final_state(bcm, hewind):
    """Collect final masses for one hewind value."""
    final_bcm = bcm.sort_values(["bin_num", "tphys"]).groupby("bin_num").last()

    summary = pd.DataFrame(index=final_bcm.index)
    summary["hewind"] = hewind
    summary["initial_he_mass_1"] = he_masses[summary.index.astype(int)]
    summary["final_mass_1"] = final_bcm["mass_1"]
    summary["mass_lost_1"] = summary["initial_he_mass_1"] - summary["final_mass_1"]
    summary["final_kstar_1"] = final_bcm["kstar_1"]
    return summary.reset_index(drop=True)


def evolve_hewind_grid(helium_grid, windflag=0, dtp=2.0):
    """Run the helium-star grid for each hewind value."""
    summaries = []

    for hewind in hewind_values:
        settings = BSEDict.copy()
        settings["windflag"] = windflag
        settings["hewind"] = hewind
        settings["random_seed"] = 1

        bpp, bcm, initC, kick_info = Evolve.evolve(
            initialbinarytable=helium_grid,
            BSEDict=settings,
            dtp=dtp,
        )
        summaries.append(summarize_final_state(bcm, hewind))

    return pd.concat(summaries, ignore_index=True)


helium_grid = make_helium_star_grid()
summaries = evolve_hewind_grid(helium_grid, windflag=0, dtp=2.0)
baseline = summaries[summaries["hewind"].eq(0.5)].set_index("initial_he_mass_1")

fig, axes = plt.subplots(1, 2, figsize=(10, 4.0), sharex=True, constrained_layout=True)

for hewind in hewind_values:
    subset = summaries[summaries["hewind"].eq(hewind)].sort_values(
        "initial_he_mass_1"
    )
    label = hewind_label(hewind)

    axes[0].plot(
        subset["initial_he_mass_1"],
        subset["final_mass_1"],
        marker="o",
        label=label,
    )

    indexed = subset.set_index("initial_he_mass_1")
    common = indexed.index.intersection(baseline.index)
    delta = indexed.loc[common, "final_mass_1"] - baseline.loc[
        common, "final_mass_1"
    ]
    axes[1].plot(common, delta, marker="o", label=label)

axes[1].axhline(0, linestyle="--", linewidth=1, color="black")
axes[0].set_xlabel("Initial He-star mass [$M_\\odot$]")
axes[0].set_ylabel("Final primary mass [$M_\\odot$]")
axes[0].set_title("Final mass")

axes[1].set_xlabel("Initial He-star mass [$M_\\odot$]")
axes[1].set_ylabel("Delta final mass [$M_\\odot$]")
axes[1].set_title("Difference from default")

for ax in axes:
    ax.legend(title="hewind", fontsize=8, title_fontsize=9, markerscale=1.2)
    ax.grid(True, alpha=0.3)

plt.savefig("hewind_flagtest.png", dpi=150)
plt.show()
