__all__ = ['ALL_COLUMNS', 'INTEGER_COLUMNS', 'BPP_COLUMNS', 'BCM_COLUMNS', 'KICK_COLUMNS', 'GROUPED_SETTINGS']

# Make this match the ordering of all_cols in bpp_array.f
ALL_COLUMNS = ['tphys', 'mass_1', 'mass_2', 'kstar_1', 'kstar_2', 'sep', 'porb',
               'ecc', 'RRLO_1', 'RRLO_2', 'evol_type', 'aj_1', 'aj_2', 'tms_1',
               'tms_2', 'massc_he_layer_1', 'massc_he_layer_2', 'massc_co_layer_1', 'massc_co_layer_2',
               'rad_1', 'rad_2', 'mass0_1',
               'mass0_2', 'lum_1', 'lum_2', 'teff_1', 'teff_2', 'radc_1',
               'radc_2', 'menv_1', 'menv_2', 'renv_1', 'renv_2', 'omega_spin_1',
               'omega_spin_2', 'B_1', 'B_2', 'bacc_1', 'bacc_2', 'tacc_1',
               'tacc_2', 'epoch_1', 'epoch_2', 'bhspin_1', 'bhspin_2',
               'deltam_1', 'deltam_2', 'SN_1', 'SN_2', 'bin_state', 'merger_type', 'metallicity']

INTEGER_COLUMNS = ["bin_state", "bin_num", "kstar_1", "kstar_2", "SN_1", "SN_2", "evol_type"]

BPP_COLUMNS = ['tphys', 'mass_1', 'mass_2', 'kstar_1', 'kstar_2',
               'sep', 'porb', 'ecc', 'RRLO_1', 'RRLO_2', 'evol_type',
               'aj_1', 'aj_2', 'tms_1', 'tms_2',
               'massc_he_layer_1', 'massc_he_layer_2', 'massc_co_layer_1', 'massc_co_layer_2', 'rad_1', 'rad_2',
               'mass0_1', 'mass0_2', 'lum_1', 'lum_2', 'teff_1', 'teff_2',
               'radc_1', 'radc_2', 'menv_1', 'menv_2', 'renv_1', 'renv_2',
               'omega_spin_1', 'omega_spin_2', 'B_1', 'B_2', 'bacc_1', 'bacc_2',
               'tacc_1', 'tacc_2', 'epoch_1', 'epoch_2',
               'bhspin_1', 'bhspin_2']

BCM_COLUMNS = ['tphys', 'kstar_1', 'mass0_1', 'mass_1', 'lum_1', 'rad_1',
               'teff_1', 'massc_he_layer_1', 'massc_co_layer_1', 'radc_1', 'menv_1', 'renv_1', 'epoch_1',
               'omega_spin_1', 'deltam_1', 'RRLO_1', 'kstar_2', 'mass0_2', 'mass_2',
               'lum_2', 'rad_2', 'teff_2', 'massc_he_layer_2', 'massc_co_layer_2', 'radc_2', 'menv_2',
               'renv_2', 'epoch_2', 'omega_spin_2', 'deltam_2', 'RRLO_2',
               'porb', 'sep', 'ecc', 'B_1', 'B_2',
               'SN_1', 'SN_2', 'bin_state', 'merger_type']

KICK_COLUMNS = ['star', 'disrupted', 'natal_kick', 'phi', 'theta', 'mean_anomaly',
                'delta_vsysx_1', 'delta_vsysy_1', 'delta_vsysz_1', 'vsys_1_total',
                'delta_vsysx_2', 'delta_vsysy_2', 'delta_vsysz_2', 'vsys_2_total',
                'theta_euler', 'phi_euler', 'psi_euler', 'randomseed', 'tphys', 'bin_num']

GROUPED_SETTINGS = {
    "windvars": [
        "neta", "bwind", "hewind", "beta", "xi", "acc2",
        "epsnov", "eddfac", "gamma", "LBV_flag",
    ],
    "cevars": ["alpha1", "lambdaf", "qcrit_array"],
    "ceflags": ["ceflag", "cekickflag", "cemergeflag", "cehestarflag", "ussn"],
    "flags": [
        "tflag", "ifflag", "wdflag", "rtmsflag", "bhflag", "remnantflag",
        "maltsev_mode", "grflag", "bhms_coll_flag", "wd_mass_lim", "aic",
        "bdecayfac", "windflag", "qcflag", "eddlimflag", "bhspinflag",
        "rejuvflag", "htpmb", "ST_cr", "ST_tide",
    ],
    "snvars": [
        "pisn", "ppi_co_shift", "ppi_extra_ml", "maltsev_fallback",
        "maltsev_pf_prob", "fryer_mass_limit", "mxns", "fryer_fmix",
        "fryer_mcrit_nsbh", "ecsn", "ecsn_mlow", "sigma", "sigmadiv",
        "bhsigmafrac", "polar_kick_angle", "natal_kick_array", "bhspinmag",
        "rembar_massloss", "kickflag", "mm_mu_ns", "mm_mu_bh",
    ],
    "points": ["pts1", "pts2", "pts3"],
    "mtvars": ["don_lim", "acc_lim", "smt_periastron_check"],
    "magvars": ["bconst", "ck"],
    "tidalvars": ["fprimc_array"],
    "mixvars": ["rejuv_fac"],
    "metvars": ["zsun"],
    "rand1": ["randomseed"]
}