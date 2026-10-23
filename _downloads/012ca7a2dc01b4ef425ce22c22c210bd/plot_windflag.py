"""
``windflag``
============

This example shows the effect of the ``windflag`` setting on massive-star
evolution. ``windflag`` selects the stellar wind mass-loss prescription used by
COSMIC. The initial metallicity and binary setup are held fixed, while only
``windflag`` is changed.

The systems are very wide binaries with low-mass companions, so the example is
focused on stellar winds rather than Roche-lobe overflow, common-envelope
evolution, or tides. The inset in the final-mass panel shows the same
comparison on log-log axes.
"""

import sys

sys.path.append("..")
import generate_default_bsedict

import matplotlib.pyplot as plt
from matplotlib.ticker import LogFormatterMathtext, LogLocator, NullFormatter
import numpy as np
import pandas as pd

from cosmic.evolve import Evolve
from cosmic.sample import InitialBinaryTable

BSEDict = generate_default_bsedict.get_default_BSE_settings(to_python=True)
plt.style.use("../_static/gallery.mplstyle")
plt.rcParams.update({
    "font.size": 9,
    "axes.titlesize": 12,
    "axes.labelsize": 10,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "legend.fontsize": 7,
    "legend.title_fontsize": 8,
    "figure.titlesize": 12,
})


windflag_labels = {
    -1: "-1 none",
    0: "0 SSE/BSE",
    1: "1 StarTrack",
    2: "2 Vink",
    3: "3 Vink+LBV",
    5: "5 0.33x",
    6: "6 Bjorklund",
    7: "7 Krticka",
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

fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.2), constrained_layout=False)
fig.subplots_adjust(left=0.08, right=0.98, bottom=0.28, top=0.82, wspace=0.38)
log_inset = axes[0].inset_axes([0.09, 0.55, 0.42, 0.36])

for offset, windflag in zip(np.linspace(-0.14, 0.14, len(windflag_values)), windflag_values):
    subset = summaries[summaries["windflag"].eq(windflag)]
    label = windflag_labels[windflag]
    log_subset = subset[
        subset["initial_mass_1"].gt(0) & subset["final_mass_1"].gt(0)
    ]
    axes[0].plot(
        log_subset["initial_mass_1"],
        log_subset["final_mass_1"],
        label=label,
    )
    log_inset.plot(
        log_subset["initial_mass_1"],
        log_subset["final_mass_1"],
        lw=1.0,
    )
    axes[1].scatter(
        subset["initial_mass_1"],
        subset["final_kstar_1"] + offset,
        s=12,
        label=label,
    )

axes[0].set_xlabel("Initial Mass [$M_\\odot$]")
axes[0].set_ylabel("Final Mass [$M_\\odot$]")
axes[0].set_title("Final mass")
axes[0].set_xscale("linear")
axes[0].set_yscale("linear")
axes[0].set_xticks([10, 20, 40, 60, 80])
axes[0].set_yticks([0, 15, 30, 45, 60])

log_positive = summaries[
    summaries["initial_mass_1"].gt(0) & summaries["final_mass_1"].gt(0)
]
log_inset.set_xscale("log", base=10)
log_inset.set_yscale("log", base=10)
log_inset.set_xlim(masses.min(), masses.max())
log_inset.set_ylim(
    log_positive["final_mass_1"].min() * 0.8,
    log_positive["final_mass_1"].max() * 1.2,
)
log_inset.xaxis.set_major_locator(LogLocator(base=10, numticks=3))
log_inset.yaxis.set_major_locator(LogLocator(base=10, numticks=3))
log_inset.xaxis.set_major_formatter(LogFormatterMathtext(base=10))
log_inset.yaxis.set_major_formatter(LogFormatterMathtext(base=10))
log_inset.xaxis.set_minor_locator(LogLocator(base=10, subs=np.arange(2, 10) * 0.1))
log_inset.yaxis.set_minor_locator(LogLocator(base=10, subs=np.arange(2, 10) * 0.1))
log_inset.xaxis.set_minor_formatter(NullFormatter())
log_inset.yaxis.set_minor_formatter(NullFormatter())
log_inset.tick_params(axis="both", which="major", labelsize=8, pad=1)
log_inset.tick_params(axis="both", which="minor", labelsize=0)
log_inset.grid(True, which="both", alpha=0.25)

axes[1].set_xlabel("Initial Mass [$M_\\odot$]")
axes[1].set_ylabel("Final kstar")
axes[1].set_title("Final kstar type")

used_kstars = sorted(summaries["final_kstar_1"].astype(int).unique())
axes[1].set_yticks(used_kstars)
axes[1].set_yticklabels([kstar_labels.get(k, str(k)) for k in used_kstars])

for ax in axes:
    ax.grid(True, alpha=0.3)

handles, labels = axes[0].get_legend_handles_labels()
fig.legend(
    handles,
    labels,
    title="windflag",
    loc="lower center",
    ncol=4,
    fontsize=7,
    title_fontsize=8,
    bbox_to_anchor=(0.5, 0.03),
)
plt.show()
