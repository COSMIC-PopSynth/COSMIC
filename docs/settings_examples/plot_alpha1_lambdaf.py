"""
alpha1 and lambdaf
==================

This example compares two common-envelope parameters using real COSMIC
separation histories.  ``alpha1`` changes the efficiency with which orbital
energy ejects the common envelope.  ``lambdaf`` changes the envelope
binding-energy lambda prescription; the negative values used here select fixed
lambda values.

Each panel holds the initial binary fixed and changes only the parameter shown
in that panel.  Open circles mark surviving binaries, x markers mark mergers,
and triangles mark disruptions at the final real COSMIC output point shown.
Tracks that overlap exactly are given tiny horizontal plotting offsets for
visibility only; the COSMIC output times and separations are unchanged.
"""

import copy
import sys
import warnings

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from cosmic.evolve import Evolve
from cosmic.sample import InitialBinaryTable

sys.path.append("..")
import generate_default_bsedict

BSEDict = generate_default_bsedict.get_default_BSE_settings(to_python=True)


EXAMPLES = [
    {
        "title": "alpha1: CE efficiency",
        "flag_key": "alpha1",
        "values": [0.2, 1.0, 5.0],
        "labels": {
            0.2: r"$\alpha$ = 0.2",
            1.0: r"$\alpha$ = 1.0",
            5.0: r"$\alpha$ = 5.0",
        },
        "initial": {
            "m1": 8.0,
            "m2": 2.8,
            "porb": 1584.893192461114,
            "tphysf": 70.0,
        },
        "plot_end_pad": 0.35,
    },
    {
        "title": "lambdaf: fixed envelope lambda",
        "flag_key": "lambdaf",
        "values": [-0.1, -0.5, -1.0],
        "labels": {
            -0.1: r"fixed $\lambda$ = 0.1",
            -0.5: r"fixed $\lambda$ = 0.5",
            -1.0: r"fixed $\lambda$ = 1.0",
        },
        "initial": {
            "m1": 13.0,
            "m2": 3.25,
            "porb": 2238.72113856834,
            "tphysf": 50.0,
        },
        "plot_end_pad": 0.5,
    },
]

LINE_WIDTHS = {
    "alpha1": {
        0.2: 3.0,
        1.0: 2.2,
        5.0: 1.4,
    },
    "lambdaf": {
        -0.1: 3.0,
        -0.5: 2.2,
        -1.0: 1.4,
    },
}

LINE_ZORDERS = {
    "alpha1": {
        0.2: 3,
        1.0: 4,
        5.0: 5,
    },
    "lambdaf": {
        -0.1: 3,
        -0.5: 4,
        -1.0: 5,
    },
}

DISPLAY_X_OFFSETS = {
    "alpha1": {
        0.2: 0.00,
        1.0: 0.15,
        5.0: 0.30,
    },
    "lambdaf": {
        -0.1: 0.00,
        -0.5: 0.15,
        -1.0: 0.30,
    },
}


def make_binary(initial):
    """Build a one-row InitialBinaryTable."""
    return InitialBinaryTable.InitialBinaries(
        m1=np.array([initial["m1"]]),
        m2=np.array([initial["m2"]]),
        porb=np.array([initial["porb"]]),
        ecc=np.zeros(1),
        tphysf=np.ones(1) * initial["tphysf"],
        kstar1=np.ones(1),
        kstar2=np.ones(1),
        metallicity=np.ones(1) * 0.014,
    )


def set_parameter(settings, key, value):
    """Set one CE parameter while preserving COSMIC's expected type."""
    if key == "alpha1":
        settings[key] = [value, value]
    else:
        settings[key] = value


def evolve_one(example, value):
    """Evolve one binary while changing only the selected CE parameter."""
    settings = copy.deepcopy(BSEDict)
    set_parameter(settings, example["flag_key"], value)

    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            message="At least one of your initial binaries is starting in Roche Lobe Overflow:*",
            category=UserWarning,
        )
        bpp, bcm, initC, kick_info = Evolve.evolve(
            initialbinarytable=make_binary(example["initial"]),
            BSEDict=settings,
            dtp=0.5,
            randomseed=1,
        )

    return bpp, bcm


def first_ce_time(bpp):
    """Return the first confirmed common-envelope event time."""
    ce_rows = bpp[bpp["evol_type"].eq(7)]
    if ce_rows.empty:
        raise RuntimeError("Selected binary does not enter common-envelope evolution.")
    return ce_rows["tphys"].iloc[0]


def terminal_time(bpp):
    """Return the first terminal-event time, or NaN if none is present."""
    terminal_rows = bpp[bpp["evol_type"].eq(6)]
    if terminal_rows.empty:
        return np.nan
    return terminal_rows["tphys"].iloc[0]


def final_status(bpp, bcm):
    """Read the final outcome from COSMIC outputs."""
    if np.isfinite(terminal_time(bpp)):
        return "merges"

    final_row = bcm.sort_values("tphys").iloc[-1]
    if final_row["bin_state"] == 0 and final_row["sep"] > 0:
        return "survives"
    if final_row["bin_state"] == 1:
        return "merges"
    return "disrupts"


def separation_track(bpp, bcm):
    """Combine sampled bcm rows with bpp event rows into a real separation track."""
    sampled_rows = bcm[["tphys", "sep"]].copy()
    sampled_rows["source_order"] = 0

    event_rows = bpp[["tphys", "sep"]].copy()
    event_rows["source_order"] = np.arange(len(event_rows)) + 1

    track = pd.concat([sampled_rows, event_rows], ignore_index=True)
    track = track[np.isfinite(track["tphys"]) & np.isfinite(track["sep"])]
    track = track[track["sep"].gt(0)]
    track = track.drop_duplicates(subset=["tphys", "sep"], keep="first")
    track = track.sort_values(["tphys", "source_order"]).reset_index(drop=True)

    stop_time = terminal_time(bpp)
    if np.isfinite(stop_time):
        track = track[track["tphys"].le(stop_time)].copy()

    assert len(track) > 0
    assert np.all(np.isfinite(track["tphys"]))
    assert np.all(np.isfinite(track["sep"]))
    assert np.all(track["sep"] > 0)
    assert np.all(np.diff(track["tphys"]) >= 0)

    return track


def mark_outcome(ax, track, status, color, x_offset=0.0, marker_scale=1.0, zorder=12):
    """Place the endpoint marker at the final real plotted point."""
    final = track.iloc[-1]
    marker_time = final["tphys"] + x_offset
    if status == "survives":
        ax.scatter(
            marker_time,
            final["sep"],
            marker="o",
            s=52 * marker_scale,
            facecolor="white",
            edgecolor=color,
            linewidths=1.8,
            zorder=zorder,
        )
    elif status == "merges":
        ax.scatter(
            marker_time,
            final["sep"],
            marker="x",
            s=62 * marker_scale,
            color=color,
            linewidths=2.1,
            zorder=zorder,
        )
    else:
        ax.scatter(
            marker_time,
            final["sep"],
            marker="^",
            s=62 * marker_scale,
            facecolor="white",
            edgecolor=color,
            linewidths=1.8,
            zorder=zorder,
        )


fig, axes = plt.subplots(
    1,
    2,
    figsize=(12, 5),
    constrained_layout=True,
)

for ax, example in zip(axes, EXAMPLES):
    colors = plt.get_cmap("tab10").colors[: len(example["values"])]
    records = []

    for value, color in zip(example["values"], colors):
        bpp, bcm = evolve_one(example, value)
        ce_time = first_ce_time(bpp)
        track = separation_track(bpp, bcm)
        status = final_status(bpp, bcm)
        records.append(
            {
                "value": value,
                "label": example["labels"][value],
                "color": color,
                "track": track,
                "status": status,
                "ce_time": ce_time,
            }
        )

    ce_times = np.array([record["ce_time"] for record in records])
    if not np.allclose(ce_times, ce_times[0]):
        raise RuntimeError("Selected runs do not share the same CE onset time.")
    ce_time = ce_times[0]

    for record in records:
        linewidth = LINE_WIDTHS[example["flag_key"]][record["value"]]
        zorder = LINE_ZORDERS[example["flag_key"]][record["value"]]
        display_offset = DISPLAY_X_OFFSETS[example["flag_key"]][record["value"]]
        marker_offset = display_offset
        marker_scale = 1.0
        marker_zorder = 12 + zorder
        if example["flag_key"] == "lambdaf" and record["status"] == "merges":
            marker_offset += {-0.1: -0.06, -0.5: 0.06}.get(record["value"], 0.0)
            marker_scale = {-0.1: 1.15, -0.5: 0.95}.get(record["value"], 1.0)
        ax.plot(
            record["track"]["tphys"] + display_offset,
            record["track"]["sep"],
            color=record["color"],
            lw=linewidth,
            linestyle="-",
            marker=None,
            label=record["label"],
            zorder=zorder,
        )
        mark_outcome(
            ax,
            record["track"],
            record["status"],
            record["color"],
            x_offset=marker_offset,
            marker_scale=marker_scale,
            zorder=marker_zorder,
        )

    combined = pd.concat([record["track"] for record in records], ignore_index=True)
    max_time = max(record["track"]["tphys"].max() for record in records)
    ax.set_xlim(0, max_time + example["plot_end_pad"])
    ax.set_ylim(combined["sep"].min() * 0.85, combined["sep"].max() * 1.35)
    ax.set_yscale("log")
    ax.axvline(ce_time, color="0.45", lw=0.9, zorder=0)
    x0, x1 = ax.get_xlim()
    ax.text(
        ce_time + 0.028 * (x1 - x0),
        0.93,
        "CE onset",
        transform=ax.get_xaxis_transform(),
        ha="left",
        va="top",
        fontsize=7,
        color="0.25",
    )
    ax.set_title(example["title"], fontsize=11)
    ax.set_xlabel("Time [Myr]")
    ax.set_ylabel(r"Orbital separation [R$_\odot$]")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(
        loc="best",
        fontsize=8,
        frameon=True,
        title=None,
        handlelength=1.8,
    )

plt.show()
