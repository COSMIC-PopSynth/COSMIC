"""
gamma
=====

This example tests ``gamma``, the angular-momentum prescription for material
lost from the system during super-Eddington Roche-lobe overflow.  The same
initial binary is evolved for each value of ``gamma``.  The selected binary
starts stable Roche-lobe overflow from the primary, loses mass from the system,
and shows a clear orbital response to the angular-momentum-loss prescription.
"""

import copy
import sys
import warnings

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

from cosmic.evolve import Evolve
from cosmic.sample import InitialBinaryTable

sys.path.append("..")
import generate_default_bsedict

BSEDict = generate_default_bsedict.get_default_BSE_settings(to_python=True)


GAMMA_VALUES = [-1, -2, -3, 0.0, 1.0, 2.0]
GAMMA_LABELS = {
    -1: "-1: primary AM",
    -2: "-2: secondary wind",
    -3: "-3: L2 disk",
    0.0: "0: no AM loss",
    1.0: "1: constant",
    2.0: "2: constant",
}

INITIAL_BINARY = {
    "m1": 3.0,
    "m2": 1.2,
    "porb": 3.0199517204020165,
    "tphysf": 500.0,
}

LINE_WIDTHS = {
    -1: 2.4,
    -2: 2.8,
    -3: 2.0,
    0.0: 2.2,
    1.0: 1.8,
    2.0: 1.6,
}


def make_binary():
    """Build the one binary used for every gamma value."""
    return InitialBinaryTable.InitialBinaries(
        m1=np.array([INITIAL_BINARY["m1"]]),
        m2=np.array([INITIAL_BINARY["m2"]]),
        porb=np.array([INITIAL_BINARY["porb"]]),
        ecc=np.zeros(1),
        tphysf=np.ones(1) * INITIAL_BINARY["tphysf"],
        kstar1=np.ones(1),
        kstar2=np.ones(1),
        metallicity=np.ones(1) * 0.014,
    )


def evolve_gamma(gamma):
    """Evolve the selected binary while changing only gamma."""
    settings = copy.deepcopy(BSEDict)
    settings["gamma"] = gamma

    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            message="At least one of your initial binaries is starting in Roche Lobe Overflow:*",
            category=UserWarning,
        )
        bpp, bcm, initC, kick_info = Evolve.evolve(
            initialbinarytable=make_binary(),
            BSEDict=settings,
            dtp=0.05,
            randomseed=1,
        )

    return bpp.sort_values("tphys").reset_index(drop=True), bcm.sort_values(
        "tphys"
    ).reset_index(drop=True)


def terminal_time(bpp):
    """Return a terminal-event time if COSMIC records one."""
    terminal_rows = bpp[bpp["evol_type"].eq(6)]
    if terminal_rows.empty:
        return np.nan
    return terminal_rows["tphys"].iloc[0]


def final_status(bpp, bcm):
    """Return the final status from COSMIC output."""
    final_row = bcm.iloc[-1]
    if final_row["bin_state"] == 0 and final_row["sep"] > 0:
        return "survives"
    if final_row["bin_state"] == 1 or np.isfinite(terminal_time(bpp)):
        return "merges"
    return "disrupts"


def first_rlof_interval(bpp):
    """Return the first RLOF start and end times."""
    starts = bpp[bpp["evol_type"].eq(3)]
    if starts.empty:
        raise RuntimeError("Selected binary does not begin Roche-lobe overflow.")

    start_time = starts["tphys"].iloc[0]
    ends = bpp[bpp["evol_type"].eq(4) & bpp["tphys"].ge(start_time)]
    if ends.empty:
        end_time = bpp["tphys"].max()
    else:
        end_time = ends["tphys"].iloc[0]

    return start_time, end_time


def rlof_diagnostics(bpp):
    """Measure non-conservative RLOF directly from the event log."""
    rlof_start, rlof_end = first_rlof_interval(bpp)
    start_row = bpp.iloc[(bpp["tphys"] - rlof_start).abs().argsort()[:1]].iloc[0]
    end_row = bpp.iloc[(bpp["tphys"] - rlof_end).abs().argsort()[:1]].iloc[0]

    donor = 1 if start_row["RRLO_1"] >= start_row["RRLO_2"] else 2
    if donor == 1:
        donor_loss = max(0.0, start_row["mass_1"] - end_row["mass_1"])
        accretor_gain = max(0.0, end_row["mass_2"] - start_row["mass_2"])
    else:
        donor_loss = max(0.0, start_row["mass_2"] - end_row["mass_2"])
        accretor_gain = max(0.0, end_row["mass_1"] - start_row["mass_1"])

    systemic_loss = max(0.0, donor_loss - accretor_gain)

    return {
        "rlof_start": rlof_start,
        "rlof_end": rlof_end,
        "donor": donor,
        "donor_loss": donor_loss,
        "accretor_gain": accretor_gain,
        "systemic_mass_loss": systemic_loss,
        "super_eddington_rlof": donor_loss > accretor_gain and systemic_loss > 1.0e-3,
    }


def evolution_track(bpp, bcm):
    """Combine sampled bcm rows and bpp event rows without inventing points."""
    sampled = bcm[["tphys", "sep", "porb"]].copy()
    sampled["source_order"] = 0
    events = bpp[["tphys", "sep", "porb"]].copy()
    events["source_order"] = np.arange(len(events)) + 1

    track = pd.concat([sampled, events], ignore_index=True)
    track = track[np.isfinite(track["tphys"]) & np.isfinite(track["sep"])]
    track = track[track["sep"].gt(0)]
    track = track.drop_duplicates(subset=["tphys", "sep", "porb"], keep="first")
    track = track.sort_values(["tphys", "source_order"]).reset_index(drop=True)

    stop_time = terminal_time(bpp)
    if np.isfinite(stop_time):
        track = track[track["tphys"].le(stop_time)].copy()

    return track


def outcome_marker(ax, track, status, color, zorder):
    """Draw only the final outcome marker."""
    final = track.iloc[-1]
    if status == "survives":
        ax.scatter(
            final["tphys"],
            final["sep"],
            marker="o",
            s=48,
            facecolor="white",
            edgecolor=color,
            linewidths=1.6,
            zorder=zorder,
        )
    elif status == "merges":
        ax.scatter(
            final["tphys"],
            final["sep"],
            marker="x",
            s=54,
            color=color,
            linewidths=2.0,
            zorder=zorder,
        )
    else:
        ax.scatter(
            final["tphys"],
            final["sep"],
            marker="^",
            s=54,
            facecolor="white",
            edgecolor=color,
            linewidths=1.6,
            zorder=zorder,
        )


def summarize(gamma, bpp, bcm, track, diagnostics):
    """Build one row of the printed summary table."""
    initial = bcm.iloc[0]
    final = bcm.iloc[-1]
    positive_sep = bcm[bcm["sep"].gt(0)]["sep"]

    return {
        "gamma": gamma,
        "super_Edd_RLOF": diagnostics["super_eddington_rlof"],
        "systemic_mass_lost": diagnostics["systemic_mass_loss"],
        "initial_sep": initial["sep"],
        "final_sep": final["sep"],
        "initial_porb": initial["porb"],
        "final_porb": final["porb"],
        "min_sep": positive_sep.min(),
        "final_kstar_1": int(final["kstar_1"]),
        "final_kstar_2": int(final["kstar_2"]),
        "status": final_status(bpp, bcm),
    }


results = []
for gamma in GAMMA_VALUES:
    bpp, bcm = evolve_gamma(gamma)
    diagnostics = rlof_diagnostics(bpp)
    track = evolution_track(bpp, bcm)
    results.append(
        {
            "gamma": gamma,
            "bpp": bpp,
            "bcm": bcm,
            "track": track,
            "diagnostics": diagnostics,
            "status": final_status(bpp, bcm),
        }
    )

if not all(item["diagnostics"]["super_eddington_rlof"] for item in results):
    raise RuntimeError("The selected binary did not show super-Eddington RLOF for all gamma values.")

if not any(item["diagnostics"]["systemic_mass_loss"] > 0 for item in results):
    raise RuntimeError("The selected binary did not lose mass from the system.")

summary = pd.DataFrame(
    [
        summarize(
            item["gamma"],
            item["bpp"],
            item["bcm"],
            item["track"],
            item["diagnostics"],
        )
        for item in results
    ]
)

rlof_start = min(item["diagnostics"]["rlof_start"] for item in results)
rlof_end = max(
    item["diagnostics"]["rlof_end"]
    for item in results
    if item["status"] == "survives"
)
x_start = max(0.0, rlof_start - 1.0)
x_end = min(
    max(item["track"]["tphys"].max() for item in results) + 1.0,
    max(rlof_end + 4.0, rlof_start + 5.0),
)

fig, ax = plt.subplots(figsize=(8.8, 5.2), constrained_layout=True)

colors = plt.get_cmap("tab10").colors[: len(GAMMA_VALUES)]

for zorder, (item, color) in enumerate(zip(results, colors), start=3):
    gamma = item["gamma"]
    track = item["track"]
    ax.plot(
        track["tphys"],
        track["sep"],
        color=color,
        lw=LINE_WIDTHS[gamma],
        linestyle="-",
        marker=None,
        label=GAMMA_LABELS[gamma],
        zorder=zorder,
    )
    outcome_marker(ax, track, item["status"], color, zorder=12 + zorder)

ax.axvline(rlof_start, color="0.35", lw=1.0)
ax.grid(True, which="both", alpha=0.25)
ax.set_xlim(x_start, x_end)
ax.set_yscale("log")
ax.set_ylabel(r"Orbital separation [$R_\odot$]")
ax.set_xlabel("Time [Myr]")
ax.set_title("gamma: angular momentum loss during super-Eddington RLOF")
gamma_legend = ax.legend(
    loc="lower right",
    fontsize=8,
    frameon=True,
    ncol=2,
    handlelength=1.8,
)
ax.add_artist(gamma_legend)
endpoint_legend = [
    Line2D(
        [0],
        [0],
        marker="o",
        color="0.25",
        markerfacecolor="white",
        lw=0,
        markersize=6,
        label="survives",
    ),
    Line2D(
        [0],
        [0],
        marker="x",
        color="0.25",
        lw=0,
        markersize=7,
        label="merges",
    ),
]
ax.legend(
    handles=endpoint_legend,
    loc="lower left",
    fontsize=8,
    frameon=True,
    handlelength=1.2,
)

y_top = ax.get_ylim()[1]
ax.text(
    rlof_start,
    y_top / 1.4,
    "RLOF starts",
    rotation=90,
    ha="right",
    va="top",
    fontsize=8,
    color="0.25",
)
plt.show()
