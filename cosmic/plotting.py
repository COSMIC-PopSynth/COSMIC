# -*- coding: utf-8 -*-
# Copyright (C) Katelyn Breivik (2017 - 2021)
#
# This file is part of cosmic.
#
# cosmic is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# cosmic is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with cosmic.  If not, see <http://www.gnu.org/licenses/>.

"""`plotting`
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches

from matplotlib.lines import Line2D
from .utils import a_from_p, calc_Roche_radius
from .evolve import Evolve

__author__ = "Jeff Andrews <andrews.jeff@gmail.com>"
__all__ = [
    "evolve_binary",
    "plot_k_type",
    "plot_radius",
    "plot_mass",
    "plot_Teff",
    "plot_Mdot",
    "plot_P_orb",
    "plot_ecc",
    "plot_HR_diagram",
    "plot_binary_evol",
    "evolve_and_plot",
]

rsun_in_au = 215.0954
day_in_year = 365.242

# Colors
primary_color = "C0"
secondary_color = "C1"


def evolve_binary(initC, t_min=None, t_max=None, BSEDict={}):
    """
    Evolves a single binary with all timesteps written
    to the bcm array for plotting

    Parameters
    ----------
    initC : `pandas.DataFrame`
        initial conditions for binary to evolve

    t_min : `float`
        starting time for plot in Myr

    t_max : `float`
        ending time for plot in Myr

    BSEDict : `Dict`
        Dictionary containing all BSE flags needed

    Returns
    -------
    bcm : `pandas.DataFrame`
        binary evolution data form BSE's bcm dict
    """

    # Disable chained warning for now which throws
    # a warning for setting dtp and tphysf
    pd.options.mode.chained_assignment = None

    # Get highest BSE temporal resolution
    BSEDict["dtp"] = 0.01

    # To deal with the limited size of the bcm array, we need to reduce
    # the time resolution if we are evolving a binary for more than 100 Myr
    if t_max is not None:
        if t_min is not None:
            if t_max - t_min > 100:
                BSEDict["dtp"] = (t_max - t_min) / 10000
        else:
            if t_max > 100:
                BSEDict["dtp"] = t_max / 10000

    # Set maximum time for evolution
    if t_max is not None:
        initC["tphysf"] = t_max

    # Call evolution scripts
    bpp, bcm, initC, kick_info = Evolve.evolve(
        initialbinarytable=initC, BSEDict=BSEDict
    )

    # In case minimum and maximum times are not supplied by user
    if t_min is None:
        # BUG HERE - BSE sets first value to weird float
        # HACKED FIX
        t_phys = np.array(bcm.tphys)
        t_phys[0] = 0.0
        bcm.at[0, "tphys"] = t_phys
        t_min = bcm.tphys.min()
    if t_max is None:
        t_max = bcm.tphys.max()

    # Select for temporal range
    bcm = bcm.loc[(bcm.tphys <= t_max) & (bcm.tphys >= t_min)]

    return bcm


def plot_k_type(ax_1, ax_2, ax_k_type_list, times, k1_out, k2_out,
                k_type_colors=None, k_type_labels=None):
    """Plots the stellar types as a function of time

    Parameters
    ----------
    ax_1 : `matplotlib.Axes`
        Axis instance for kstar 1

    ax_2 : `matplotlib.Axes`
        Axis instance for kstar 2

    ax_k_type_list : `list`
        List of ktypes for the legend of the ktype bar plots

    times : `pandas.Series`
        Series of times at each binary evolution time step in Myr

    k1_out : `pandas.Series`
        Series of kstar 1 type at each binary evolution time step in Myr

    k2_out : `pandas.Series`
        Series of kstar 2 type at each binary evolution time step in Myr

    k_type_colors : `list`
        List of colors for each ktype

    k_type_labels : `list`
        List of labels for each ktype

    Returns
    -------
    Null
    """

    k_type_colors = [
        "plum",
        "sandybrown",
        "lightseagreen",
        "moccasin",
        "chartreuse",
        "deepskyblue",
        "gold",
        "rosybrown",
        "m",
        "darkgreen",
        "grey",
        "sienna",
        "palevioletred",
        "navy",
        "tan",
        "black",
    ] if k_type_colors is None else k_type_colors
    k_type = [
        "MS conv",
        "MS",
        "HG",
        "GB",
        "CHeB",
        "EAGB",
        "TPAGB",
        "HeMS",
        "HeHG",
        "HeGB",
        "HeWD",
        "COWD",
        "ONeWD",
        "NS",
        "BH",
        "no remnant",
    ] if k_type_labels is None else k_type_labels

    # k-type plots
    for a in [ax_1, ax_2]:

        a.axis("off")

        for j, k in enumerate(np.unique(k1_out)):
            a.fill_between(times[k1_out == k], 0.7, 0.9, color=k_type_colors[k])
        for j, k in enumerate(np.unique(k2_out)):
            a.fill_between(times[k2_out == k], 0.4, 0.6, color=k_type_colors[k])

        a.set_title("Stellar Type")

        a.set_xlim(np.min(times), np.max(times))

    # Add legend
    ax_k_type_list.axis("off")
    k_type_all = np.unique(np.array([k1_out, k2_out]))
    patches = [
        mpatches.Patch(color=k_type_colors[k], label=k_type[k]) for k in k_type_all
    ]
    leg = ax_k_type_list.legend(handles=patches, mode="expand", ncol=len(k_type_all))
    leg.get_frame().set_alpha(0.0)


def plot_radius(
    ax, times, R1_out, R2_out, M1_out, M2_out, P_orb_out, ecc_out, sys_obs, legend=True
):
    """Plots the stellar radii as a function of time

    Parameters
    ----------
    ax : `matplotlib.Axes`
        Axis instance for the plot

    times : `pandas.Series`
        Series of times at each binary evolution time step in Myr

    R1_out : `pandas.Series`
        Series of primary radius at each binary evolution time step in Rsun

    R2_out : `pandas.Series`
        Series of secondary radius at each binary evolution time step in Rsun

    M1_out : `pandas.Series`
        Series of primary mass at each binary evolution time step in Msun

    M2_out : `pandas.Series`
        Series of secondary mass at each binary evolution time step in Msun

    P_orb_out : `pandas.Series`
        Series of orbital periods at each binary evolution time step in days

    ecc_out : `pandas.Series`
        Series of eccentricities at each binary evolution time step

    sys_obs : `Dict`
        Dictionary containing keys for binary parameters with values to plot
        as vertical lines for each stellar component

    Returns
    -------
    Null
    """
    # Radius
    ax.plot(times, R1_out, color=primary_color)
    ax.plot(times, R2_out, color=secondary_color)

    # Plot Roche radii - at periastron
    A_out = a_from_p(M1_out, M2_out, P_orb_out)  # in Rsun
    R1_Roche = calc_Roche_radius(M1_out, M2_out, A_out * (1.0 - ecc_out))  # in Rsun
    R2_Roche = calc_Roche_radius(M2_out, M1_out, A_out * (1.0 - ecc_out))  # in Rsun
    ax.plot(times, R1_Roche, color=primary_color, linestyle="--")
    ax.plot(times, R2_Roche, color=secondary_color, linestyle="--")

    for key, value in sys_obs.items():
        if key == "R1":
            ax.axhline(value, color=primary_color, linestyle="dashed")
        if key == "R2":
            ax.axhline(value, color=secondary_color, linestyle="dashed")

    ax.set_yscale("log")
    ax.set_ylim(0.5, ax.get_ylim()[1])

    ax.set_ylabel(r"Radius ($R_{\odot}$)")

    ax.set_xlim(np.min(times), np.max(times))

    if legend:
        custom_lines = [
            Line2D([0], [0], color="k", linestyle="solid"),
            Line2D([0], [0], color="k", linestyle="dashed"),
        ]
        ax.legend(custom_lines, ["Stellar Radius", "Roche Radius"], frameon=False)


def plot_mass(ax, times, M1_out, M2_out, sys_obs):
    """Plots the stellar mass as a function of time

    Parameters
    ----------
    ax : `matplotlib.Axes`
        Axis instance for the plot

    times : `pandas.Series`
        Series of times at each binary evolution time step in Myr

    M1_out : `pandas.Series`
        Series of primary mass at each binary evolution time step in Msun

    M2_out : `pandas.Series`
        Series of secondary mass at each binary evolution time step in Msun

    sys_obs : `Dict`
        Dictionary containing keys for binary parameters with values to plot
        as vertical lines for each stellar component

    Returns
    -------
    Null
    """

    # Mass
    ax.plot(times, M1_out, color=primary_color)
    ax.plot(times, M2_out, color=secondary_color)

    for key, value in sys_obs.items():
        if key == "M1":
            ax.axhline(value, color=primary_color, linestyle="dashed")
        if key == "M2":
            ax.axhline(value, color=secondary_color, linestyle="dashed")

    ax.set_ylabel(r"Mass ($M_{\odot}$)")
    ax.set_xlim(np.min(times), np.max(times))


def plot_Teff(ax, times, Teff1_out, Teff2_out, sys_obs):
    """Plots the stellar effectvie temperature as a function of time

    Parameters
    ----------
    ax : `matplotlib.Axes`
        Axis instance for the plot

    times : `pandas.Series`
        Series of times at each binary evolution time step in Myr

    Teff1_out : `pandas.Series`
        Series of primary effective temperature at each binary evolution time step in K

    Teff2_out : `pandas.Series`
        Series of secondary effective temperature at each binary evolution time step in K

    sys_obs : `Dict`
        Dictionary containing keys for binary parameters with values to plot
        as vertical lines for each stellar component

    Returns
    -------
    Null
    """

    # Effective temperature
    ax.plot(times, Teff1_out, color=primary_color)
    ax.plot(times, Teff2_out, color=secondary_color)

    for key, value in sys_obs.items():
        if key == "T1":
            ax.axhline(value, color=primary_color, linestyle="dashed")
        if key == "T2":
            ax.axhline(value, color=secondary_color, linestyle="dashed")

    ax.set_yscale("log")
    ax.set_ylabel(r"T$_{\rm eff}$ (K)")
    ax.set_xlim(np.min(times), np.max(times))
    ax.set_xlabel("Time (Myr)")


def plot_Mdot(ax, times, Mdot1_out, Mdot2_out, legend=True):
    """Plots the stellar effectvie temperature as a function of time

    Parameters
    ----------
    ax : `matplotlib.Axes`
        Axis instance for the plot

    times : `pandas.Series`
        Series of times at each binary evolution time step in Myr

    Mdot1_out : `pandas.Series`
        Series of primary mass transfer rate at each binary evolution time step in Msun/yr

    Mdot2_out : `pandas.Series`
        Series of secondary mass transfer rate at each binary evolution time step in Msun/yr

    sys_obs : `Dict`
        Dictionary containing keys for binary parameters with values to plot
        as vertical lines for each stellar component

    Returns
    -------
    Null
    """

    # Mass accretion rate
    ax.plot(times, np.log10(np.clip(Mdot1_out, 1.0e-16, None)), color=primary_color)
    ax.plot(times, np.log10(np.clip(Mdot2_out, 1.0e-16, None)), color=secondary_color)
    ax.set_ylabel(r"Mass Accretion Rate ($M_{\odot}$ yr$^{-1}$)")

    ax.set_ylim(-14, ax.get_ylim()[1])
    ax.set_xlim(np.min(times), np.max(times))

    if legend:
        custom_lines = [
            Line2D([0], [0], color="C0", linestyle="solid"),
            Line2D([0], [0], color="C1", linestyle="solid"),
        ]
        ax.legend(custom_lines, ["Primary", "Secondary"], loc=2, frameon=False)


def plot_P_orb(ax, times, P_orb_out, t_max, sys_obs):
    """Plots the orbital period as a function of time

    Parameters
    ----------
    ax : `matplotlib.Axes`
        Axis instance for the plot

    times : `pandas.Series`
        Series of times at each binary evolution time step in Myr

    P_orb_out : `pandas.Series`
        Series of orbital periods at each binary evolution time step in days

    t_max : `float`
        ending time for plot in Myr

    sys_obs : `Dict`
        Dictionary containing keys for binary parameters with values to plot
        as vertical lines for each stellar component

    Returns
    -------
    Null
    """

    ax.plot(times, P_orb_out, color="k")

    for key, value in sys_obs.items():
        if key == "P_orb":
            ax.axhline(value, color="k", linestyle="dashed")

    ax.set_xlim(np.min(times), t_max)

    ax.set_ylabel(r"P$_{\rm orb}$ (days)")
    ax.set_yscale("log")


def plot_ecc(ax, times, ecc_out, t_max, sys_obs):
    """Plots the eccentricity as a function of time

    Parameters
    ----------
    ax : `matplotlib.Axes`
        Axis instance for the plot

    times : `pandas.Series`
        Series of times at each binary evolution time step in Myr

    ecc_out : `pandas.Series`
        Series of eccentricities at each binary evolution time step

    t_max : `float`
        ending time for plot in Myr

    sys_obs : `Dict`
        Dictionary containing keys for binary parameters with values to plot
        as vertical lines for each stellar component

    Returns
    -------
    Null
    """

    ax.plot(times, ecc_out, color="k")

    for key, value in sys_obs.items():
        if key == "ecc":
            ax.axhline(value, color="k", linestyle="dashed")

    ax.set_xlim(np.min(times), t_max)

    ax.set_ylabel("Eccentricity")
    ax.set_xlabel("Time (Myr)")


def plot_HR_diagram(ax, L1_out, L2_out, Teff1_out, Teff2_out):
    """Plots the eccentricity as a function of time

    Parameters
    ----------
    ax : `matplotlib.Axes`
        Axis instance for the plot

    times : `pandas.Series`
        Series of times at each binary evolution time step in Myr

    L1_out : `pandas.Series`
        Series of primary luminosities at each binary evolution time step in Lsun

    L2_out : `pandas.Series`
        Series of secondary luminosities at each binary evolution time step in Lsun

    Teff1_out : `pandas.Series`
        Series of primary effective temperature at each binary evolution time step in K

    Teff2_out : `pandas.Series`
        Series of secondary effective temperature at each binary evolution time step in K

    sys_obs : `Dict`
        Dictionary containing keys for binary parameters with values to plot
        as vertical lines for each stellar component

    Returns
    -------
    Null
    """

    ax.plot(Teff1_out, L1_out, color=primary_color)
    ax.plot(Teff2_out, L2_out, color=secondary_color)

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(40000, 2000)

    ax.set_xlabel(r"log T$_{\rm eff}$ (K)")
    ax.set_ylabel("log L (erg/s)")

    ax.set_xticks([2000, 5000, 10000, 20000, 40000])


def plot_binary_evol(bcm, sys_obs={}, ktype_kwargs={}, t_min=None, t_max=None):
    """Plots the full set of plots and kstar bar plots

    Parameters
    ----------
    bcm : `pandas.DataFrame`
        binary evolution data form BSE's bcm dict

    sys_obs : `Dict`
        Dictionary containing keys for binary parameters with values to plot
        as vertical lines for each stellar component

    Returns
    -------
    fig : `matplotlib.figure`
        Figure containing all the plots!
    """

    fig = plt.figure(figsize=(10, 10), layout="tight")

    gs = fig.add_gridspec(2, 1, height_ratios=[1.3, 6], hspace=0.0)
    gs_top = gs[0].subgridspec(2, 2, width_ratios=[1, 1], height_ratios=[1, 1.3], hspace=0.15)
    gs_bottom = gs[1].subgridspec(3, 2, width_ratios=[1, 1], height_ratios=[2, 2, 2], hspace=0.15)

    if t_min is not None:
        bcm = bcm.loc[bcm.tphys >= t_min]
    if t_max is not None:
        bcm = bcm.loc[bcm.tphys <= t_max]

    # k-type panels
    plot_k_type(
        fig.add_subplot(gs_top[0, 0]),
        fig.add_subplot(gs_top[0, 1]),
        fig.add_subplot(gs_top[1, :]),
        bcm.tphys,
        bcm.kstar_1.astype(int),
        bcm.kstar_2.astype(int),
        **ktype_kwargs
    )

    # Radius panel
    plot_radius(
        fig.add_subplot(gs_bottom[0, 0]),
        bcm.tphys,
        bcm.rad_1,
        bcm.rad_2,
        bcm.mass_1,
        bcm.mass_2,
        bcm.porb,
        bcm.ecc,
        sys_obs,
    )

    # Mass panel
    plot_mass(fig.add_subplot(gs_bottom[1, 0]), bcm.tphys, bcm.mass_1, bcm.mass_2, sys_obs)

    # Teff panel
    plot_Teff(fig.add_subplot(gs_bottom[2, 0]), bcm.tphys, bcm.teff_1, bcm.teff_2, sys_obs)

    # Mass accretion rate panel
    plot_Mdot(fig.add_subplot(gs_bottom[0, 1]), bcm.tphys, bcm.deltam_1, bcm.deltam_2)

    # Orbital period panel
    plot_P_orb(
        fig.add_subplot(gs_bottom[1, 1]),
        bcm.loc[bcm.kstar_2 < 15].tphys,
        bcm.loc[bcm.kstar_2 < 15].porb,
        np.max(bcm.tphys),
        sys_obs,
    )

    # Plot eccentricity panel
    plot_ecc(
        fig.add_subplot(gs_bottom[2, 1]),
        bcm.loc[bcm.kstar_2 < 15].tphys,
        bcm.loc[bcm.kstar_2 < 15].ecc,
        np.max(bcm.tphys),
        sys_obs,
    )

    # Plot HR diagram
    # idx_1 = np.where(k1_out < 10)[0]
    # idx_2 = np.where(k2_out < 10)[0]
    # plot_HR_diagram(ax[7], L1_out[k1_out<10], L2_out[k2_out<10], Teff1_out[k1_out<10], Teff2_out[k2_out<10])

    return fig


def evolve_and_plot(initC, t_min=None, t_max=None, BSEDict=None, sys_obs={}):
    """
    Evolve and plot binaries as a function of time

    Parameters
    ----------
    initC : `pandas.DataFrame`
        initial conditions for binary to evolve

    t_min : `float or list`
        starting time for plot in Myr

    t_max : `float or list`
        ending time for plot in Myr

    BSEDict : `Dict`
        Dictionary containing all BSE flags needed

    sys_obs : `Dict`
        Dictionary containing keys for binary parameters with values to plot
        as vertical lines for each stellar component

    Returns
    -------
    all_figs : `list`
        list of all figures created
    """

    # Throw an error if user is plotting more than 10 binaries
    if len(initC) > 10:
        raise ValueError(
            "You have asked to plot more than 10 separate binaries. This could cause problems..."
        )

    # Iterate over all binaries in initC
    all_figs = []
    for i in range(len(initC)):

        # Check if t_min and t_max are lists
        if isinstance(t_min, list):
            t_min_tmp = t_min[i]
        else:
            t_min_tmp = t_min
        if isinstance(t_max, list):
            t_max_tmp = t_max[i]
        else:
            t_max_tmp = t_max

        # Evolve binary
        bcm = evolve_binary(
            initC.iloc[i: i + 1], t_min=t_min_tmp, t_max=t_max_tmp, BSEDict=BSEDict
        )

        # Plot binary
        fig = plot_binary_evol(bcm, sys_obs=sys_obs)

        # Add to list of figs
        all_figs.append(fig)

    # Return list of figs
    return all_figs
