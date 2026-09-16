"""
``bwind``
=========

This example tests ``bwind``, the binary-enhanced wind mass-loss parameter from
Hurley et al. (2000), Equation 12.  The same initial binary is evolved for each
value of ``bwind``.  The selected system remains detached, but the primary
expands close enough to its Roche lobe that the binary-enhanced wind term
changes the mass loss before ordinary Roche-lobe overflow or common-envelope
evolution can dominate.
"""

import copy
import sys
import warnings

import matplotlib.pyplot as plt
import numpy as np

from cosmic.evolve import Evolve
from cosmic.sample import InitialBinaryTable

sys.path.append("..")
import generate_default_bsedict

try:
    plt.style.use("../_static/gallery.mplstyle")
except OSError:
    pass

BSEDict = generate_default_bsedict.get_default_BSE_settings(to_python=True)


BWIND_VALUES = [0.0, 1000.0, 10000.0]
BWIND_LABELS = {
    0.0: "0: none",
    1000.0: "1,000",
    10000.0: "10,000",
}

SEARCH_INITIALS = {
    "m1": 1.5,
    "m2": 1.0,
    "periods": np.array([1800.0, 2100.0, 2346.2288481422624, 2700.0, 3100.0]),
}

FINAL_TPHYS = 14000.0
SEARCH_TPHYS = 4000.0
PLOT_DTP = 0.1
METALLICITY = 0.014
RLOF_EVOL_TYPE = 3
CE_EVOL_TYPE = 7
TERMINAL_EVOL_TYPE = 6
CLOSE_RRLO_THRESHOLD = 0.25


def make_binary(m1, m2, porb, tphysf):
    """Build one initial binary."""
    return InitialBinaryTable.InitialBinaries(
        m1=np.array([m1]),
        m2=np.array([m2]),
        porb=np.array([porb]),
        ecc=np.zeros(1),
        tphysf=np.ones(1) * tphysf,
        kstar1=np.ones(1),
        kstar2=np.ones(1),
        metallicity=np.ones(1) * METALLICITY,
    )


def evolve_bwind(m1, m2, porb, bwind, tphysf=FINAL_TPHYS, dtp=1.0):
    """Evolve one binary while changing only bwind."""
    settings = copy.deepcopy(BSEDict)
    settings["bwind"] = bwind

    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            message="At least one of your initial binaries is starting in Roche Lobe Overflow:*",
            category=UserWarning,
        )
        bpp, bcm, initC, kick_info = Evolve.evolve(
            initialbinarytable=make_binary(m1, m2, porb, tphysf),
            BSEDict=settings,
            dtp=dtp,
            randomseed=1,
        )

    return (
        bpp.sort_values("tphys").reset_index(drop=True),
        bcm.sort_values("tphys").reset_index(drop=True),
    )


def interaction_flags(bpp, bcm):
    """Return whether RLOF, CE, or merger-like termination occurred."""
    final = bcm.iloc[-1]
    return {
        "rlof": bool(bpp["evol_type"].eq(RLOF_EVOL_TYPE).any()),
        "ce": bool(bpp["evol_type"].eq(CE_EVOL_TYPE).any()),
        "merger": bool(
            bpp["evol_type"].eq(TERMINAL_EVOL_TYPE).any()
            or final["bin_state"] != 0
            or final["sep"] <= 0
        ),
    }


def close_phase_rows(bcm):
    """Rows where the primary is close enough to its Roche lobe for bwind."""
    close = bcm[bcm["RRLO_1"].ge(CLOSE_RRLO_THRESHOLD)].copy()
    if close.empty:
        peak_index = bcm["RRLO_1"].idxmax()
        close = bcm.loc[[peak_index]].copy()
    return close


def quick_metrics(bpp, bcm):
    """Metrics used by the small search."""
    close = close_phase_rows(bcm)
    initial = bcm.iloc[0]
    final = bcm.iloc[-1]
    flags = interaction_flags(bpp, bcm)
    return {
        **flags,
        "max_rrlo": float(bcm["RRLO_1"].max()),
        "mass_lost": float(initial["mass_1"] - final["mass_1"]),
        "close_mass_lost": float(close["mass_1"].iloc[0] - close["mass_1"].iloc[-1]),
        "max_kstar_1": int(bcm["kstar_1"].max()),
    }


def select_binary():
    """Run a small period search and return the first clean bwind-sensitive case."""
    m1 = SEARCH_INITIALS["m1"]
    m2 = SEARCH_INITIALS["m2"]

    for porb in SEARCH_INITIALS["periods"]:
        metrics = []
        for bwind in BWIND_VALUES:
            bpp, bcm = evolve_bwind(m1, m2, porb, bwind, SEARCH_TPHYS, dtp=2.0)
            metrics.append(quick_metrics(bpp, bcm))

        if any(row["rlof"] or row["ce"] or row["merger"] for row in metrics):
            continue
        if max(row["max_kstar_1"] for row in metrics) < 3:
            continue
        if max(row["max_rrlo"] for row in metrics) < 0.49:
            continue
        if max(row["max_rrlo"] for row in metrics) >= 1.0:
            continue
        if max(row["mass_lost"] for row in metrics) - min(row["mass_lost"] for row in metrics) < 0.01:
            continue

        return {"m1": m1, "m2": m2, "porb": porb}

    raise RuntimeError("No clean detached bwind-sensitive binary found in the search grid.")


def choose_time_window(results):
    """Zoom around the giant phase where at least one model approaches RL filling."""
    active_times = []
    for bcm in results.values():
        active = bcm[bcm["RRLO_1"].ge(0.08)]
        if not active.empty:
            active_times.extend([active["tphys"].min(), active["tphys"].max()])

    if not active_times:
        return 0.0, max(bcm["tphys"].max() for bcm in results.values())

    start = max(0.0, min(active_times) - 150.0)
    end = min(max(bcm["tphys"].max() for bcm in results.values()), max(active_times) + 250.0)
    return start, end


selected = select_binary()

evolved = {}
for bwind in BWIND_VALUES:
    bpp, bcm = evolve_bwind(
        selected["m1"], selected["m2"], selected["porb"], bwind, FINAL_TPHYS, dtp=PLOT_DTP
    )
    evolved[bwind] = {"bpp": bpp, "bcm": bcm}

primary_mass_lost = [
    row["bcm"]["mass_1"].iloc[0] - row["bcm"]["mass_1"].iloc[-1]
    for row in evolved.values()
]
if max(primary_mass_lost) - min(primary_mass_lost) < 0.01:
    raise RuntimeError("The selected binary does not show a measurable bwind mass-loss effect.")

fig, axes = plt.subplots(3, 1, figsize=(9.0, 9.2), sharex=True)
colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
line_widths = [3.4, 2.6, 1.8]
# The three RRLO_1 peaks land within ~1 Myr of each other, so a single wide,
# top-zorder curve there would visually swallow the shorter peaks beneath it.
# Draw that panel with thinner, equal-width lines in decreasing-peak-height
# order (tallest first/back, shortest last/front) so all three show up as
# nested spikes instead.
RRLO_LINE_WIDTH = 2.0

for i, bwind in enumerate(BWIND_VALUES):
    bcm = evolved[bwind]["bcm"]
    color = colors[i % len(colors)]
    lw = line_widths[i % len(line_widths)]
    label = BWIND_LABELS[bwind]

    axes[0].plot(
        bcm["tphys"],
        bcm["mass_1"],
        color=color,
        lw=lw,
        label=label,
        zorder=10 - i,
    )
    axes[1].plot(
        bcm["tphys"],
        bcm["RRLO_1"],
        color=color,
        lw=RRLO_LINE_WIDTH,
        label=label,
        zorder=10 + i,
    )
    axes[2].plot(
        bcm["tphys"],
        bcm["sep"],
        color=color,
        lw=lw,
        label=label,
        zorder=10 - i,
    )

x_min, _ = choose_time_window({value: row["bcm"] for value, row in evolved.items()})
x_max = 3150.0

axes[0].set_ylabel(r"Primary mass [$M_\odot$]", fontsize=12)
axes[1].set_ylabel(r"Roche-lobe filling factor $R_1/R_{\mathrm{L},1}$", fontsize=12)
axes[2].set_ylabel(r"Separation [$R_\odot$]", fontsize=12)
axes[2].set_xlabel("Time [Myr]", fontsize=12)
max_rrlo = max(row["bcm"]["RRLO_1"].max() for row in evolved.values())
axes[1].set_ylim(0.0, max(0.75, max_rrlo * 1.12))

for ax in axes:
    ax.set_xlim(x_min, x_max)
    ax.grid(alpha=0.25)
    ax.tick_params(labelsize=10)

axes[0].legend(
    title="bwind",
    loc="best",
    fontsize=9,
    title_fontsize=9,
    frameon=True,
)

axes[0].set_title("Binary-enhanced winds near Roche-lobe filling", fontsize=14)

fig.tight_layout()
plt.show()
