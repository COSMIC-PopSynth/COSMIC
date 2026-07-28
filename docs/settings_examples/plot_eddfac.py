""" 
eddfac 
======

This example tests ``eddfac``, the Eddington-limit for accretion onto a compact object like a NS, BH, or WD: the accretor can accrete at up to ``eddfac`` times the Eddington rate beforre the excess transferred mass is lost from the system. ``efffac=0`` disables cmopact-object accretion entirely; the default is ``eddfact=1``, and this uses ``eddfac=10``. 

A pre-formed 1.4 Msun NS accretes from a range of donor masses and orbital
periods (Z=0.02), swept over ``eddfac`` in {0, 1, 3, 10, 30, 100}. The left
panel shows the final accretor mass for the most accretion-sensitive
donor/period combinations; the right panel shows the mean and max mass
accreted across all surviving systems, as a function of ``eddfac``.
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from cosmic.sample.initialbinarytable import InitialBinaryTable
from cosmic.evolve import Evolve

BSEDict = {
    'pts1': 0.001, 'pts2': 0.01, 'pts3': 0.02, 'zsun': 0.014, 'windflag': 3,
    'eddlimflag': 0, 'neta': 0.5, 'bwind': 0.0, 'hewind': 0.5, 'beta': 0.125,
    'xi': 0.5, 'acc2': 1.5, 'LBV_flag': 1, 'alpha1': [1.0, 1.0], 'lambdaf': 0.0,
    'ceflag': 1, 'cekickflag': 2, 'cemergeflag': 1, 'cehestarflag': 0, 'qcflag': 5,
    'qcrit_array': [0.0]*16, 'kickflag': 5, 'sigma': 265.0, 'bhflag': 1,
    'bhsigmafrac': 1.0, 'sigmadiv': -20.0, 'ecsn': 2.25, 'ecsn_mlow': 1.6,
    'aic': 1, 'ussn': 1, 'polar_kick_angle': 90.0,
    'natal_kick_array': [[-100.0, -100.0, -100.0, -100.0, 0.0], [-100.0, -100.0, -100.0, -100.0, 0.0]],
    'mm_mu_ns': 400.0, 'mm_mu_bh': 200.0, 'remnantflag': 4, 'fryer_mass_limit': 0,
    'mxns': 3.0, 'fryer_fmix': 1.0, 'fryer_mcrit_nsbh': 5.75, 'rembar_massloss': 0.5,
    'wd_mass_lim': 1, 'maltsev_mode': 0, 'maltsev_fallback': 0.5, 'maltsev_pf_prob': 0.1,
    'pisn': -2, 'ppi_co_shift': 0.0, 'ppi_extra_ml': 0.0, 'bhspinflag': 0,
    'bhspinmag': 0.0, 'grflag': 1, 'eddfac': 10, 'gamma': -2, 'don_lim': -1,
    'acc_lim': [-1, -1], 'smt_periastron_check': 0, 'tflag': 1, 'ST_tide': 1,
    'fprimc_array': [2.0/21.0]*16, 'ifflag': 1, 'wdflag': 1, 'epsnov': 0.001,
    'bdecayfac': 1, 'bconst': 3000, 'ck': 1000, 'rejuv_fac': 1.0, 'rejuvflag': 0,
    'bhms_coll_flag': 0, 'htpmb': 1, 'ST_cr': 1, 'rtmsflag': 0,
}

outdir = os.path.expanduser('~/projects/')
Z = 0.02
m2g = [1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0]
pbg = [10, 30, 100, 300, 1000]
M2, PB = np.meshgrid(m2g, pbg)
m2arr = M2.ravel().astype(float)
pbarr = PB.ravel().astype(float)
n = len(m2arr)
eddvals = [0.0, 1.0, 3.0, 10.0, 30.0, 100.0]

def run(edd):
    BSEDict['eddfac'] = edd
    ibt = InitialBinaryTable.InitialBinaries(
        m1=[1.4]*n, m2=list(m2arr), porb=list(pbarr), ecc=[0.0]*n,
        tphysf=[13700.0]*n, kstar1=[13]*n, kstar2=[1]*n, metallicity=[Z]*n,
    )
    bpp, bcm, initC, kick_info = Evolve.evolve(initialbinarytable=ibt, BSEDict=BSEDict)
    fin = bcm.reset_index().groupby('bin_num').last()
    return fin['mass_1'], fin['kstar_1']

mass_tab, kstar_tab = {}, {}
for edd in eddvals:
    m, k = run(edd)
    mass_tab[edd], kstar_tab[edd] = m, k
    d = m.values - 1.4
    print(f'eddfac={edd:6.1f}: {int((d>0.01).sum()):2d} accreting, '
          f'max gain {d.max():.3f} Msun, {int((k.values==14).sum()):2d} -> BH', flush=True)

Mass = pd.DataFrame(mass_tab)
Kst = pd.DataFrame(kstar_tab)
survive = (Kst.isin([13, 14])).all(axis=1)
Mass = Mass[survive]
spread = (Mass.max(axis=1) - Mass.min(axis=1))
top = spread.sort_values(ascending=False).head(6).index

# use the actual eddfac values on the x-axis (symlog so eddfac=0 can still be shown)
x = np.array(eddvals)

fig, ax = plt.subplots(1, 2, figsize=(13, 5.2))
for c, b in zip(plt.cm.viridis(np.linspace(0.1, 0.9, len(top))), top):
    ax[0].plot(x, Mass.loc[b, eddvals].values, marker='o', lw=1.8, color=c,
               label=fr'$M_2$={m2arr[b]:.1f}, $P$={pbarr[b]:.0f}d')
ax[0].axhline(1.4, ls=':', color='0.5', lw=1.0)
ax[0].set_xscale('symlog', linthresh=1)
ax[0].set_xticks(eddvals)
ax[0].set_xticklabels([f'{e:g}' for e in eddvals])
ax[0].set_xlim(-0.5, 100)
ax[0].set_ylim(1.38, Mass.loc[top, eddvals].values.max() + 0.05)
ax[0].set_xlabel('eddfac')
ax[0].set_ylabel(r'Final accretor mass [$M_\odot$]')
ax[0].set_title('eddfac: final accretor mass')
ax[0].legend(frameon=False, fontsize=11)
ax[0].grid(alpha=0.25, which='both')

dmean = (Mass.values - 1.4).mean(axis=0)
dmax = (Mass.values - 1.4).max(axis=0)
ax[1].plot(x, dmean, marker='o', lw=1.8, color='C0', label='mean over surviving systems')
ax[1].plot(x, dmax, marker='s', lw=1.8, color='C3', label='max over surviving systems')
ax[1].set_xscale('symlog', linthresh=1)
ax[1].set_xticks(eddvals)
ax[1].set_xticklabels([f'{e:g}' for e in eddvals])
ax[1].set_xlim(-0.5, 100)
ax[1].set_xlabel('eddfac')
ax[1].set_ylabel(r'Mass accreted $\Delta M$ [$M_\odot$]')
ax[1].set_title('eddfac: accreted mass')
ax[1].legend(frameon=False, fontsize=11)
ax[1].grid(alpha=0.25, which='both')

plt.tight_layout()
plt.savefig(os.path.join(outdir, 'eddfac_test.png'), dpi=160)
plt.show()
print('saved:', os.path.join(outdir, 'eddfac_test.png'), flush=True)
