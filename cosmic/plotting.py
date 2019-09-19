# -*- coding: utf-8 -*-
# Copyright (C) Katelyn Breivik (2017 - 2019)
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

'''`plotting`
'''

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
from .utils import a_from_p, p_from_a, calc_Roche_radius
from .evolve import Evolve


__author__ = 'Jeff Andrews <andrews.jeff@gmail.com>'
__all__ = ['evolve_binary', 'plot_k_type', 'plot_radius',
           'plot_mass','plot_Teff','plot_Mdot','plot_P_orb','plot_ecc',
           'plot_HR_diagram','plot_binary_evol','evolve_and_plot']


rsun_in_au = 215.0954
day_in_year = 365.242


# Colors
primary_color = 'C0'
secondary_color = 'C1'


def evolve_binary(initC, t_min=None, t_max=None, BSEDict={}):
    """
    Docstring to be added.
    """
    # Disable chained warning for now which throws
    # a warning for setting dtp and tphysf
    pd.options.mode.chained_assignment = None 

    # Get highest BSE temporal resolution
    initC['dtp'] = 0.01

    # Set maximum time for evolution
    if not t_max is None: initC['tphysf'] = t_max

    # Call evolution scripts
    bpp, bcm, initC = Evolve.evolve(initialbinarytable=initC, BSEDict=BSEDict)

    # In case minimum and maximum times are not supplied by user
    if t_min is None:
        # BUG HERE - BSE sets first value to weird float
        # HACKED FIX
        t_phys = np.array(bcm.tphys)
        t_phys[0] = 0.0
        bcm.at[0, 'tphys'] = t_phys
        t_min = bcm.tphys.min()
    if t_max is None: t_max = bcm.tphys.max()

    # Select for temporal range
    bcm = bcm.loc[(bcm.tphys<=t_max) & (bcm.tphys>=t_min)]

    return bcm


def plot_k_type(ax_1, ax_2, ax_k_type_list, times, k1_out, k2_out):


    k_type_colors = ['plum', 'sandybrown', 'lightseagreen', 'moccasin', 'chartreuse',
                     'deepskyblue', 'gold', 'rosybrown', 'm', 'darkgreen', 'grey',
                     'sienna', 'palevioletred', 'navy', 'tan', 'black']
    k_type = ['MS conv', 'MS', 'HG', 'GB', 'CHeB', 'EAGB', 'TPAGB', 'HeMS', 'HeHG',
              'HeGB', 'HeWD', 'COWD', 'ONeWD', 'NS', 'BH', 'no remnant']

    # k-type plots
    for a in [ax_1, ax_2]:

        a.axis('off')

        for j,k in enumerate(np.unique(k1_out)):
            a.fill_between(times[k1_out == k], 0.7, 0.9, color=k_type_colors[k])
        for j,k in enumerate(np.unique(k2_out)):
            a.fill_between(times[k2_out == k], 0.4, 0.6, color=k_type_colors[k])

        a.set_title('Stellar Type')

        a.set_xlim(np.min(times), np.max(times))



    # Add legend
    ax_k_type_list.axis('off')
    k_type_all = np.unique(np.array([k1_out, k2_out]))
    patches = [ mpatches.Patch(color=k_type_colors[k], label=k_type[k]) for k in k_type_all ]
    leg = ax_k_type_list.legend(handles=patches, mode='expand', ncol=len(k_type_all))
    leg.get_frame().set_alpha(0.0)



def plot_radius(ax, times, R1_out, R2_out, M1_out, M2_out, P_orb_out, ecc_out, sys_obs, legend=True):


    # Radius
    ax.plot(times, R1_out, color=primary_color)
    ax.plot(times, R2_out, color=secondary_color)

    # Plot Roche radii - at periastron
    A_out = a_from_p(M1_out, M2_out, P_orb_out) # in Rsun
    R1_Roche = calc_Roche_radius(M1_out, M2_out, A_out*(1.0-ecc_out)) # in Rsun
    R2_Roche = calc_Roche_radius(M2_out, M1_out, A_out*(1.0-ecc_out)) # in Rsun
    ax.plot(times, R1_Roche, color=primary_color, linestyle='--')
    ax.plot(times, R2_Roche, color=secondary_color, linestyle='--')


    for key, value in sys_obs.items():
        if key == 'R1': ax.axhline(value, color=primary_color, linestyle='dashed')
        if key == 'R2': ax.axhline(value, color=secondary_color, linestyle='dashed')


    ax.set_yscale('log')
    ax.set_ylim(0.5, ax.get_ylim()[1])

    ax.set_ylabel(r'Radius ($R_{\odot}$)')

    ax.set_xlim(np.min(times), np.max(times))

    if legend:
        custom_lines = [Line2D([0], [0], color='k', linestyle='solid'),
                        Line2D([0], [0], color='k', linestyle='dashed')]
        ax.legend(custom_lines, ['Stellar Radius', 'Roche Radius'], frameon=False)


def plot_mass(ax, times, M1_out, M2_out, sys_obs):

    # Mass
    ax.plot(times, M1_out, color=primary_color)
    ax.plot(times, M2_out, color=secondary_color)

    for key, value in sys_obs.items():
        if key == 'M1': ax.axhline(value, color=primary_color, linestyle='dashed')
        if key == 'M2': ax.axhline(value, color=secondary_color, linestyle='dashed')

    ax.set_ylabel(r'Mass ($M_{\odot}$)')
    ax.set_xlim(np.min(times), np.max(times))



def plot_Teff(ax, times, Teff1_out, Teff2_out, sys_obs):



    # Effective temperature
    ax.plot(times, Teff1_out, color=primary_color)
    ax.plot(times, Teff2_out, color=secondary_color)

    for key, value in sys_obs.items():
        if key == 'T1': ax.axhline(value, color=primary_color, linestyle='dashed')
        if key == 'T2': ax.axhline(value, color=secondary_color, linestyle='dashed')

    ax.set_yscale('log')
    ax.set_ylabel(r'T$_{\rm eff}$ (K)')
    ax.set_xlim(np.min(times), np.max(times))
    ax.set_xlabel("Time (Myr)")


def plot_Mdot(ax, times, Mdot1_out, Mdot2_out, legend=True):

    # Mass accretion rate
    ax.plot(times, np.log10(np.clip(Mdot1_out, 1.0e-16, None)), color=primary_color)
    ax.plot(times, np.log10(np.clip(Mdot2_out, 1.0e-16, None)), color=secondary_color)
    ax.set_ylabel(r'Mass Accretion Rate ($M_{\odot}$ yr$^{-1}$)')

    ax.set_ylim(-14, ax.get_ylim()[1])
    ax.set_xlim(np.min(times), np.max(times))

    if legend:
        custom_lines = [Line2D([0], [0], color='C0', linestyle='solid'),
                        Line2D([0], [0], color='C1', linestyle='solid')]
        ax.legend(custom_lines, ['Primary', 'Secondary'], loc=2, frameon=False)



def plot_P_orb(ax, times, P_orb_out, t_max, sys_obs):

    ax.plot(times, P_orb_out*day_in_year, color='k')

    for key, value in sys_obs.items():
        if key == 'P_orb': ax.axhline(value, color='k', linestyle='dashed')

    ax.set_xlim(np.min(times), t_max)

    ax.set_ylabel(r'P$_{\rm orb}$ (days)')
    ax.set_yscale('log')


def plot_ecc(ax, times, ecc_out, t_max, sys_obs):

    ax.plot(times, ecc_out, color='k')

    for key, value in sys_obs.items():
        if key == 'ecc': ax.axhline(value, color='k', linestyle='dashed')


    ax.set_xlim(np.min(times), t_max)

    ax.set_ylabel('Eccentricity')
    ax.set_xlabel("Time (Myr)")


def plot_HR_diagram(ax, L1_out, L2_out, Teff1_out, Teff2_out):

    ax.plot(Teff1_out, L1_out, color=primary_color)
    ax.plot(Teff2_out, L2_out, color=secondary_color)


    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(40000, 2000)

    ax.set_xlabel(r"log T$_{\rm eff}$ (K)")
    ax.set_ylabel("log L (erg/s)")

    ax.set_xticks([2000, 5000, 10000, 20000, 40000])



def plot_binary_evol(bcm, file_out=None, sys_obs={}):


    fig, ax = plt.subplots(5, 1, figsize=(10,10))


    gs = gridspec.GridSpec(4, 2,
                           width_ratios=[2,2],
                           height_ratios=[1,2,2,2]
                           )

    ax = np.array([plt.subplot(x) for x in gs])


    # k-type panels
    ax_k_type_list = fig.add_axes([0.08, 0.75, 0.9, 0.1])
    plot_k_type(ax[0], ax[1], ax_k_type_list, bcm.tphys, bcm.kstar_1.astype(int), bcm.kstar_2.astype(int))

    # Radius panel
    plot_radius(ax[2], bcm.tphys, bcm.rad_1, bcm.rad_2, bcm.mass_1, bcm.mass_2, bcm.porb, bcm.ecc, sys_obs)

    # Mass panel
    plot_mass(ax[4], bcm.tphys, bcm.mass_1, bcm.mass_2, sys_obs)

    # Teff panel
    plot_Teff(ax[6], bcm.tphys, bcm.teff_1, bcm.teff_2, sys_obs)

    # Mass accretion rate panel
    plot_Mdot(ax[3], bcm.tphys, bcm.deltam_1, bcm.deltam_2)

    # Orbital period panel
    plot_P_orb(ax[5], bcm.loc[bcm.kstar_2<15].tphys, bcm.loc[bcm.kstar_2<15].porb, np.max(bcm.tphys), sys_obs)

    # Plot eccentricity panel
    plot_ecc(ax[7], bcm.loc[bcm.kstar_2<15].tphys, bcm.loc[bcm.kstar_2<15].ecc, np.max(bcm.tphys), sys_obs)


    # Plot HR diagram
    # idx_1 = np.where(k1_out < 10)[0]
    # idx_2 = np.where(k2_out < 10)[0]
    # plot_HR_diagram(ax[7], L1_out[k1_out<10], L2_out[k2_out<10], Teff1_out[k1_out<10], Teff2_out[k2_out<10])

    # make the labels look nice
    gs.tight_layout(fig)

    return fig




def evolve_and_plot(initC, t_min=None, t_max=None, BSEDict=None, file_out=None, sys_obs={}):
    """
    Evolve and plot binaries. Detailed docstring to be added.
    """

    # Throw an error if user is plotting more than 10 binaries
    if len(initC) > 10:
        raise ValueError("You have asked to plot more than 10 separate binaries. This could cause problems...")

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
        bcm = evolve_binary(initC.iloc[i:i+1], t_min=t_min_tmp, t_max=t_max_tmp, BSEDict=BSEDict)

        # Plot binary
        fig = plot_binary_evol(bcm, file_out=file_out, sys_obs=sys_obs)

        # Add to list of figs
        all_figs.append(fig)

    # Return list of figs
    return all_figs
