import numpy as np
from cosmic import Evolve
from cosmic.sample.initialbinarytable import InitialBinaryTable

# Create initial parameters for a single star system
# We'll create a binary where one star has negligible mass to effectively simulate a single star
initial_binaries = InitialBinaryTable.InitialBinaries(
    m1=np.array([1.0]),    # 1 solar mass primary
    m2=np.array([0.001]),  # negligible mass secondary
    porb=np.array([1e6]),  # very wide orbit to minimize interaction
    ecc=np.array([0.0]),   # circular orbit
    kstar1=np.array([1]),  # main sequence star primary
    kstar2=np.array([1]),  # main sequence star secondary
    metallicity=np.array([0.02]),  # solar metallicity
    tphysf=np.array([13700.0]),   # evolve for Hubble time
)

# Set up BSE Dict with htpmb=2
bse_dict = {
    'xi': 1.0,
    'bhflag': 1,
    'neta': 0.5,
    'windflag': 3,
    'wdflag': 1,
    'alpha1': 1.0,
    'pts1': 0.001,
    'pts3': 0.02,
    'pts2': 0.01,
    'epsnov': 0.001,
    'hewind': 1.0,
    'ck': 1000,
    'bwind': 0.0,
    'lambdaf': 1.0,
    'mxns': 3.0,
    'beta': -1.0,
    'tflag': 1,
    'acc2': 1.5,
    'remnantflag': 3,
    'ceflag': 0,
    'eddfac': 1.0,
    'ifflag': 0,
    'bconst': 3000,  # Changed from -3000 to 3000 as bconst must be positive
    'sigma': 265.0,
    'gamma': -2.0,
    'pisn': 45.0,
    'natal_kick_array': [[-100.0,-100.0,-100.0,-100.0,0.0], [-100.0,-100.0,-100.0,-100.0,0.0]],
    'bhsigmafrac': 1.0,
    'polar_kick_angle': 90,
    'qcrit_array': [0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],
    'cekickflag': 2,
    'cehestarflag': 0,
    'cemergeflag': 0,
    'ecsn': 2.25,
    'ecsn_mlow': 1.6,
    'aic': 1,
    'ussn': 0,
    'sigmadiv': -20.0,
    'qcflag': 2,
    'eddlimflag': 0,
    'fprimc_array': [2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0],
    'bhspinflag': 0,
    'bhspinmag': 0.0,
    'rejuv_fac': 1.0,
    'rejuvflag': 0,
    'htpmb': 2,  # Using our new htpmb=2 flag
    'ST_cr': 1,
    'ST_tide': 0,
    'bdecayfac': 1,
    'rembar_massloss': 0.5,
    'kickflag': 1,  # Changed from 0 to 1 as kickflag must be between 1-5
    'zsun': 0.014,
    'grflag': 1,
    'bhms_coll_flag': 0,
    'don_lim': -1,
    'acc_lim': -1,
    'rtmsflag': 0,      # Added: Controls if stars can grow beyond their MS radius limit
    'wd_mass_lim': 1.0, # Added: Mass limit for WD formation
}

# Create an Evolve object
bpp, bcm, initC, kick_info = Evolve.evolve(
    initialbinarytable=initial_binaries, 
    BSEDict=bse_dict,
    nproc=1,
)

print("\nEvolution completed!")
print("\nInitial conditions:")
print(initC)
print("\nDetailed Binary Evolution Parameters (BCM):")
print("Time (Myr) Mass_1 (Msun) Luminosity_1 (Lsun) Radius_1 (Rsun) Temperature_1 (K)")
print("-" * 80)
for i in range(len(bcm)):
    print(f"{bcm['tphys'].iloc[i]:.2f} {bcm['mass_1'].iloc[i]:.3f} {bcm['lum_1'].iloc[i]:.3e} {bcm['rad_1'].iloc[i]:.3f} {bcm['teff_1'].iloc[i]:.0f}")

print("\nKey Evolutionary Events (BPP):")
print("Time (Myr) Event_Type Mass_1 (Msun) Star_Type_1")
print("-" * 50)
for i in range(len(bpp)):
    print(f"{bpp['tphys'].iloc[i]:.2f} {bpp['evol_type'].iloc[i]:10d} {bpp['mass_1'].iloc[i]:.3f} {bpp['kstar_1'].iloc[i]:d}")

print("\nKick Information:")
print(kick_info)