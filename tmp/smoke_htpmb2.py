import numpy as np
from cosmic.evolve import Evolve
from cosmic.sample.initialbinarytable import InitialBinaryTable

# minimal single binary
m1 = [30.0]
m2 = [20.0]
porb = [10.0]
ecc = [0.1]
tphysf = [10.0]
kstar1 = [1]
kstar2 = [1]
metallicity = [0.014]

init = InitialBinaryTable.InitialBinaries(m1, m2, porb, ecc, tphysf, kstar1, kstar2, metallicity)
# set required BSE params minimally
init['neta'] = 0.5
init['bwind'] = 0.0
init['hewind'] = 0.5
init['alpha1'] = 1.0
init['lambdaf'] = 0.0
init['ceflag'] = 1
init['tflag'] = 0
init['ifflag'] = 0
init['wdflag'] = 0
init['pisn'] = -2
init['rtmsflag'] = 0
init['bhflag'] = 1
init['remnantflag'] = 4
init['grflag'] = 0
init['bhms_coll_flag'] = 0
init['wd_mass_lim'] = 1
init['cekickflag'] = 2
init['cemergeflag'] = 1
init['cehestarflag'] = 0
init['mxns'] = 3.0
init['pts1'] = 0.001
init['pts2'] = 0.01
init['pts3'] = 0.02
# set htpmb to 2
init['htpmb'] = 2
# fill other required defaults
init['randomseed'] = np.array([42])
init['bin_num'] = np.array([0])
init['natal_kick_array'] = [[[-100.0,-100.0,-100.0,-100.0,0],[-100.0,-100.0,-100.0,-100.0,0]]]
init['qcrit_array'] = [[0.0]*16]
init['fprimc_array'] = [[0.0]*16]

# run evolve: pass BSEDict instead of adding all BSE params as columns
BSEDict = {
	'neta': 0.5,
	'bwind': 0.0,
	'hewind': 0.5,
	'alpha1': 1.0,
	'lambdaf': 0.0,
	'ceflag': 1,
	'tflag': 0,
	'ifflag': 0,
	'wdflag': 0,
	'pisn': -2,
	'rtmsflag': 0,
	'bhflag': 1,
	'remnantflag': 4,
	'grflag': 0,
	'bhms_coll_flag': 0,
	'wd_mass_lim': 1,
	'cekickflag': 2,
	'cemergeflag': 1,
	'cehestarflag': 0,
	'mxns': 3.0,
	'pts1': 0.001,
	'pts2': 0.01,
	'pts3': 0.02,
	'htpmb': 2,
}

# ensure initial table contains required BSE columns that evolve() indexes
from cosmic.evolve import INITIAL_CONDITIONS_PASS_COLUMNS
defaults = {
	'ecsn': 2.25, 'ecsn_mlow': 1.6, 'aic': 1, 'ussn': 1,
	'sigma': 265.0, 'sigmadiv': -20.0, 'bhsigmafrac': 1.0, 'polar_kick_angle': 90.0,
	'beta': -1, 'xi': 0.5, 'acc2': 1.5, 'epsnov': 0.0,
	'eddfac': 1.0, 'gamma': 1.0, 'don_lim': -1, 'acc_lim': -1,
	'bdecayfac': 1, 'bconst': 1.0, 'ck': 1.0, 'windflag': 3,
	'qcflag': 1, 'eddlimflag': 0, 'bhspinflag': 0, 'bhspinmag': 0.0,
	'rejuv_fac': 1.0, 'rejuvflag': 0, 'ST_cr': 1, 'ST_tide': 1,
	'rembar_massloss': 0.5, 'zsun': 0.014, 'kickflag': 1,
}

for col in INITIAL_CONDITIONS_PASS_COLUMNS:
	if col not in init.columns:
		if col in BSEDict:
			init[col] = np.array([BSEDict[col]])
		elif col in defaults:
			init[col] = np.array([defaults[col]])
		elif col == 'natal_kick_array':
			init[col] = [[-100.0,-100.0,-100.0,-100.0,0],[-100.0,-100.0,-100.0,-100.0,0]]
		elif col == 'qcrit_array':
			init[col] = [0.0]*16
		elif col == 'fprimc_array':
			init[col] = [0.0]*16
		else:
			# generic zero
			init[col] = np.zeros(len(init))

	# ensure all columns required for saving are present as well
	from cosmic.evolve import INITIAL_BINARY_TABLE_SAVE_COLUMNS
	for col in INITIAL_BINARY_TABLE_SAVE_COLUMNS:
		if col not in init.columns:
			if col.startswith('qcrit_') or col.startswith('fprimc_'):
				init[col] = 0.0
			elif col.startswith('natal_kick_') or col.startswith('phi_') or col.startswith('theta_') or col.startswith('mean_anomaly_') or col.startswith('randomseed_'):
				init[col] = 0.0
			elif col == 'binfrac':
				init[col] = np.ones(len(init))
			else:
				init[col] = np.zeros(len(init))

	bpp, bcm, init_out, kick_info = Evolve.evolve(init, nproc=1, BSEDict=BSEDict)

print('Evolve returned: bpp shape', None if bpp is None else bpp.shape)
print('Evolve returned: bcm shape', None if bcm is None else bcm.shape)
print('Success')
print('bpp shape:', None if bpp is None else bpp.shape)
print('bcm shape:', None if bcm is None else bcm.shape)
print('Success')
