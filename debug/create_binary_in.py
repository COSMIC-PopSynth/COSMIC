import pandas as pd
import argparse

BSE_settings = {'xi': 1.0, 'bhflag': 1, 'neta': 0.5, 'windflag': 3, 'wdflag': 1, 'alpha1': 1.0,
                'pts1': 0.001, 'pts3': 0.02, 'pts2': 0.01, 'epsnov': 0.001, 'hewind': 0.5,
                'ck': 1000, 'bwind': 0.0, 'lambdaf': 0.0, 'mxns': 3.0, 'beta': -1.0, 'tflag': 1,
                'acc2': 1.5, 'grflag': 1, 'remnantflag': 4, 'ceflag': 0, 'eddfac': 1.0,
                'ifflag': 0, 'bconst': 3000, 'sigma': 265.0, 'gamma': -2.0, 'pisn': 45.0,
                'natal_kick_array': [[-100.0, -100.0, -100.0, -100.0, 0.0],
                                     [-100.0, -100.0, -100.0, -100.0, 0.0]], 'bhsigmafrac': 1.0,
                'polar_kick_angle': 90, 'qcrit_array': [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                                        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                'cekickflag': 2, 'cehestarflag': 0, 'cemergeflag': 0, 'ecsn': 2.25,
                'ecsn_mlow': 1.6, 'aic': 1, 'ussn': 0, 'sigmadiv': -20.0, 'qcflag': 5,
                'eddlimflag': 0, 'fprimc_array': [2.0/21.0, 2.0/21.0, 2.0/21.0, 2.0/21.0,
                                                  2.0/21.0, 2.0/21.0, 2.0/21.0, 2.0/21.0,
                                                  2.0/21.0, 2.0/21.0, 2.0/21.0, 2.0/21.0,
                                                  2.0/21.0, 2.0/21.0, 2.0/21.0, 2.0/21.0],
                'bhspinflag': 0, 'bhspinmag': 0.0, 'rejuv_fac': 1.0, 'rejuvflag': 0, 'htpmb': 1,
                'ST_cr': 1, 'ST_tide': 1, 'bdecayfac': 1, 'rembar_massloss': 0.5, 'kickflag' : 1,
                'zsun': 0.014, 'bhms_coll_flag': 0, 'don_lim': -1, 'acc_lim': -1, 'binfrac': 0.5,
                'rtmsflag': 0, 'wd_mass_lim': 1, 'idum': 100}


def create_binary_in(mass0, tphysf, tb, kstar, Z, ecc, BSE_settings):
    """Create a binary.in file based on the given parameters

    This follows the format in cosmic/src/test_bse.f and changes there would need to be reflected here.
    """
    with open('binary.in', 'w') as f:
        f.write(f'{mass0[0]} {mass0[1]} {tphysf} {tb} {kstar[0]} {kstar[1]} {Z} {ecc}\n')

        lines = [
            ['neta', 'bwind', 'hewind', 'alpha1', 'lambdaf', 'windflag', 'rtmsflag'],
            ['ceflag', 'tflag', 'ifflag', 'wdflag', 'bhflag', 'remnantflag', 'mxns', 'idum'],
            ['pts1', 'pts2', 'pts3'],
            ['sigma', 'beta', 'xi', 'acc2', 'epsnov', 'eddfac', 'gamma', 'kickflag'],
            ['pisn', 'cekickflag', 'cehestarflag', 'grflag', 'bhms_coll_flag'],
            ['wd_mass_lim', 'ecsn', 'ecsn_mlow', 'aic', 'ussn', 'sigmadiv', 'bhsigmafrac'],
            ['don_lim', 'acc_lim', 'bdecayfac', 'bconst', 'ck', 'qcflag', 'eddlimflag'],
            ['bhspinflag', 'bhspinmag', 'rejuv_fac', 'rejuvflag', 'htpmb', 'ST_cr'],
            ['ST_tide', 'rembar_massloss', 'zsun']
        ]

        for line in lines:
            f.write(' '.join([str(BSE_settings[key]) for key in line]) + '\n')


def convert_initC_row_to_binary_in(initC_file, bin_num):
    """Convert a row from an initC file to a binary.in file

    Parameters
    ----------
    initC_file : `str`
        Path to the initC file
    bin_num : `int`
        The binary number to convert
    """
    # get binary from initC
    initC = pd.read_hdf(initC_file, key="initC")
    r = initC.loc[bin_num]

    # update BSE settings with those in the binary
    BSE_settings['idum'] = r['randomseed'].astype(int)
    for key in BSE_settings:
        if key in r:
            BSE_settings[key] = r[key].astype(type(BSE_settings[key]))

    # create binary.in file
    create_binary_in([r['mass_1'], r['mass_2']], r['tphysf'], r['porb'],
                     [r['kstar_1'].astype(int), r['kstar_2'].astype(int)],
                     r['metallicity'], r['ecc'], BSE_settings)


def main():
    parser = argparse.ArgumentParser(description='Convert a row from an initC file to a binary.in file')
    parser.add_argument('-f', '--initC_file', type=str, help='Path to the initC file')
    parser.add_argument('-b', '--bin_num', type=int, help='The binary number to convert')
    args = parser.parse_args()

    convert_initC_row_to_binary_in(args.initC_file, args.bin_num)


if __name__ == '__main__':
    main()
