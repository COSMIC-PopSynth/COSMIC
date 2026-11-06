"""Convert a row of an initC table into the inputs needed to debug COSMIC with test_bse.f

This writes ``binary.in`` (and ``metisse_tracks.bin`` when the binary uses METISSE) in the output
directory, which ``debug/test`` reads when it runs. The values are passed to the fortran in the same
way as ``cosmic.evolve`` does, so the debug executable should follow the same evolution.

``binary.in`` uses fortran namelists (see ``src/cosmic/src/test_bse.f``), which are matched by name,
so if a setting is added to ``GROUPED_SETTINGS`` it also needs to be added to the namelist in
``test_bse.f`` (otherwise the fortran will stop with an error when reading ``binary.in``).

For METISSE, the tracks are read with the same python code that ``Evolve.evolve`` uses, but rather
than passing them to the fortran we record what would have been passed and write it to
``metisse_tracks.bin`` (see ``src/cosmic/src/test_metisse_tracks.f90`` for the reader).

Example
-------
    python create_binary_in.py -f ../bad_initC.csv -b 1479
"""
import argparse
import math
import os
from types import SimpleNamespace
from unittest import mock
import warnings

import numpy as np
import pandas as pd

from cosmic import evolve
from cosmic.consts import GROUPED_SETTINGS

# settings that are arrays in the fortran, and the initC columns that they are built from
ARRAY_SETTINGS = {
    'natal_kick_array': (evolve.FLATTENED_NATAL_KICK_COLUMNS, (2, len(evolve.NATAL_KICK_COLUMNS))),
    'qcrit_array': (evolve.QCRIT_COLUMNS, (16,)),
    'fprimc_array': (evolve.FPRIMC_COLUMNS, (16,)),
    'alpha1': (evolve.ALPHA_COLUMNS, (2,)),
    'acc_lim': (evolve.ACCLIM_COLUMNS, (2,)),
}

# the arguments that cosmic.evolve passes to evolv2: fortran name -> initC column(s)
BINARY_IN = {
    'kstar': ['kstar_1', 'kstar_2'],
    'mass0': ['mass0_1', 'mass0_2'],
    'mass': ['mass_1', 'mass_2'],
    'tb': 'porb',
    'ecc': 'ecc',
    'z': 'metallicity',
    'tphysf': 'tphysf',
    'dtp': 'dtp',
    'rad': ['rad_1', 'rad_2'],
    'lum': ['lum_1', 'lum_2'],
    'massc': ['massc_1', 'massc_2'],
    'radc': ['radc_1', 'radc_2'],
    'menv': ['menv_1', 'menv_2'],
    'renv': ['renv_1', 'renv_2'],
    'ospin': ['omega_spin_1', 'omega_spin_2'],
    'B_0': ['B_1', 'B_2'],
    'bacc': ['bacc_1', 'bacc_2'],
    'tacc': ['tacc_1', 'tacc_2'],
    'epoch': ['epoch_1', 'epoch_2'],
    'tms': ['tms_1', 'tms_2'],
    'bhspin': ['bhspin_1', 'bhspin_2'],
    'tphys': 'tphys',
}

# arguments of c_m_interface.set_format_controls_h/he, in the order of the fortran subroutine
H_FORMAT_INTS = ['read_eep', 'prems', 'zams', 'iams', 'tams', 'bgb', 'cheign', 'cheburn', 'ta_cheb',
                 'tpagb', 'ccburn', 'postagb', 'initeep', 'finaleep', 'fixtrack', 'loweep', 'higheep']
HE_FORMAT_INTS = ['read_eep', 'bgb', 'cheburn', 'ta_cheb', 'tpagb', 'ccburn', 'postagb',
                  'initeep', 'finaleep', 'fixtrack', 'loweep', 'higheep']
FORMAT_STRS = ['age_col', 'mass_col', 'logl_col', 'logt_col', 'logr_col', 'he_mass_col', 'co_mass_col',
               'he_radius_col', 'co_radius_col', 'mass_env_col', 'radius_env_col', 'logtc_col',
               'he4_col', 'c12_col', 'o16_col']
F2PY_STRLEN = 256


def fortran_value(v):
    """Format a python value as a fortran namelist value"""
    if isinstance(v, (bool, np.bool_)):
        return 'T' if v else 'F'
    if isinstance(v, str):
        return "'" + v.replace("'", "''") + "'"
    v = float(v)
    if math.isnan(v):
        return 'NaN'
    if math.isinf(v):
        return 'Inf' if v > 0 else '-Inf'
    # write integers without a decimal point so they can be read into INTEGER variables too
    if v.is_integer() and abs(v) < 2**31:
        return str(int(v))
    return repr(v)


def namelist(name, values):
    """Format a dict of {fortran name: value(s)} as a fortran namelist group"""
    lines = [f'&{name}']
    for key, val in values.items():
        val = np.asarray(val, dtype=object) if not isinstance(val, str) else val
        if isinstance(val, str) or val.ndim == 0:
            lines.append(f'  {key} = {fortran_value(val if isinstance(val, str) else val.item())}')
        elif val.ndim == 1:
            lines.append(f'  {key}(1:{len(val)}) = ' + ', '.join(fortran_value(v) for v in val))
        else:
            # one line per row, i.e. natal_kick_array(i,1:5) holds the kick for the i-th SN
            for i, row in enumerate(val):
                lines.append(f'  {key}({i + 1},1:{len(row)}) = ' + ', '.join(fortran_value(v) for v in row))
    lines.append('/')
    return '\n'.join(lines) + '\n'


def is_blank(v):
    """Whether an initC value is empty (empty strings come back as NaN from a csv)"""
    return v is None or (isinstance(v, float) and math.isnan(v)) or (isinstance(v, str) and v.strip() == '')


def get_binary_in(r):
    """Get the initial state of the binary, as passed to evolv2 by cosmic.evolve"""
    r = r.copy()
    # these mirror the defaults used by cosmic.evolve when the columns aren't in the initC
    for i in [1, 2]:
        if f'mass0_{i}' not in r:
            r[f'mass0_{i}'] = r[f'mass_{i}']
    if 'dtp' not in r:
        r['dtp'] = r['tphysf']

    values, missing = {}, []
    for key, cols in BINARY_IN.items():
        cols_list = cols if isinstance(cols, list) else [cols]
        if not all(c in r for c in cols_list):
            missing.extend(c for c in cols_list if c not in r)
            continue
        values[key] = [r[c] for c in cols] if isinstance(cols, list) else r[cols]
    if missing:
        warnings.warn(f"initC is missing {missing}, these will be set to zero in test_bse.f")
    return values


def get_bse_settings(r):
    """Get the BSE settings, as set in the fortran common blocks by cosmic.evolve"""
    values, missing = {}, []
    for settings in GROUPED_SETTINGS.values():
        for setting in settings:
            # special case: we call randomseed idum1 in the fortran
            name = 'idum1' if setting == 'randomseed' else setting.lower()
            if setting in ARRAY_SETTINGS:
                cols, shape = ARRAY_SETTINGS[setting]
                if not all(c in r for c in cols):
                    missing.append(setting)
                    continue
                values[name] = np.array([r[c] for c in cols], dtype=float).reshape(shape)
            elif setting in r:
                values[name] = r[setting]
            else:
                missing.append(setting)
    if missing:
        warnings.warn(f"initC is missing settings {missing}, these will be set to zero in test_bse.f")
    return values


def get_stellar_engine(r):
    """Get the stellar engine settings, as set in the fortran common blocks by cosmic.evolve"""
    engine = r.get('stellar_engine', 'sse')
    engine = 'sse' if is_blank(engine) else engine
    if engine not in ['sse', 'metisse']:
        raise ValueError("Use either 'sse' or 'metisse' as stellar engine")

    values = {
        'using_sse': int(engine == 'sse'),
        'using_metisse': int(engine == 'metisse'),
        'path_to_tracks': '',
        'path_to_he_tracks': '',
        'metisse_verbose': False,
    }
    if 'z_accuracy_limit' in r:
        values['z_match_limit'] = r['z_accuracy_limit']
    if engine == 'metisse':
        for col in ['path_to_tracks', 'path_to_he_tracks']:
            values[col] = '' if is_blank(r.get(col)) else str(r[col])
    return values


class _CMInterfaceRecorder:
    """Stands in for ``_evolvebin.c_m_interface``, recording what python passes to the fortran"""
    def __init__(self):
        self.format_controls = {}
        self.tracks = {}

    def set_format_controls_h(self, **kwargs):
        self.format_controls[False] = kwargs

    def set_format_controls_he(self, **kwargs):
        self.format_controls[True] = kwargs

    def set_tracks_from_python(self, *args):
        self.tracks[bool(args[-1])] = args


def record_metisse_tracks(stellar_engine, z):
    """Run cosmic.evolve.read_tracks_for_METISSE, recording what it would pass to the fortran"""
    recorder = _CMInterfaceRecorder()
    with mock.patch.object(evolve, '_evolvebin', SimpleNamespace(c_m_interface=recorder)):
        m_min, m_max = evolve.read_tracks_for_METISSE(path_to_tracks=stellar_engine['path_to_tracks'],
                                                      IBT_Z=z,
                                                      z_accuracy_limit=stellar_engine['z_match_limit'],
                                                      is_he=False)
        if stellar_engine['path_to_he_tracks'] != '':
            evolve.read_tracks_for_METISSE(path_to_tracks=stellar_engine['path_to_he_tracks'],
                                           IBT_Z=z,
                                           z_accuracy_limit=stellar_engine['z_match_limit'],
                                           is_he=True)
    return recorder, m_min, m_max


def write_metisse_tracks(recorder, path):
    """Write the recorded METISSE inputs to a stream file for test_metisse_tracks.f90"""
    has_he = True in recorder.tracks
    with open(path, 'wb') as f:
        f.write(b'METTRK01')
        f.write(np.int32(has_he).tobytes())
        for is_he in ([False, True] if has_he else [False]):
            # format controls (logicals are written as 0/1)
            fmt = recorder.format_controls[is_he]
            int_keys = HE_FORMAT_INTS if is_he else H_FORMAT_INTS
            if set(fmt) != set(int_keys) | set(FORMAT_STRS):
                raise ValueError("The arguments passed to set_format_controls_h/he have changed, "
                                 "update create_binary_in.py and test_metisse_tracks.f90 to match")
            f.write(np.int32(len(int_keys)).tobytes())
            f.write(np.array([int(fmt[k]) for k in int_keys], dtype=np.int32).tobytes())
            f.write(b''.join(str(fmt[k]).encode('ascii')[:F2PY_STRLEN].ljust(F2PY_STRLEN)
                             for k in FORMAT_STRS))

            # tracks, exactly as they would be passed to set_tracks_from_python
            (filenames, initial_mass, initial_Y, initial_Z, Fe_div_H, alpha_div_Fe, v_div_vcrit,
             ntrack_arr, neep_arr, ncol_arr, eep_data, tr_data, col_names, _) = recorder.tracks[is_he]
            f.write(np.array([len(filenames), tr_data.shape[0], eep_data.shape[0], tr_data.shape[1]],
                             dtype=np.int32).tobytes())
            f.write(np.asarray(filenames, dtype=f'S{F2PY_STRLEN}').tobytes())
            for arr in [initial_mass, initial_Y, initial_Z, Fe_div_H, alpha_div_Fe, v_div_vcrit]:
                f.write(np.asarray(arr, dtype=np.float64).tobytes())
            for arr in [ntrack_arr, neep_arr, ncol_arr]:
                f.write(np.asarray(arr, dtype=np.int32).tobytes())
            f.write(np.asarray(eep_data, dtype=np.int32).tobytes(order='F'))
            f.write(np.asarray(tr_data, dtype=np.float64).tobytes(order='F'))
            f.write(np.asarray(col_names).tobytes(order='F'))


def read_initC_row(initC_file, bin_num=None):
    """Read a single binary from an initC file (hdf5 with key initC, or csv)"""
    if initC_file.endswith('.csv'):
        initC = pd.read_csv(initC_file)
    else:
        initC = pd.read_hdf(initC_file, key="initC")

    if bin_num is None:
        if len(initC) != 1:
            raise ValueError(f"initC has {len(initC)} binaries, choose one with --bin_num")
        return initC.iloc[0]

    rows = initC[initC['bin_num'] == bin_num] if 'bin_num' in initC else initC.loc[[bin_num]]
    if len(rows) != 1:
        raise ValueError(f"Found {len(rows)} binaries with bin_num={bin_num} in {initC_file}")
    return rows.iloc[0]


def convert_initC_row_to_binary_in(initC_file, bin_num=None, out_dir='.'):
    """Convert a row from an initC file to a binary.in file (and metisse_tracks.bin if needed)

    Parameters
    ----------
    initC_file : `str`
        Path to the initC file
    bin_num : `int`
        The binary number to convert (can be omitted if the initC has a single binary)
    out_dir : `str`
        Directory in which to write the files (where debug/test is run)
    """
    r = read_initC_row(initC_file, bin_num)

    binary_in = get_binary_in(r)
    bse_settings = get_bse_settings(r)
    stellar_engine = get_stellar_engine(r)

    with open(os.path.join(out_dir, 'binary.in'), 'w') as f:
        f.write(namelist('binary_in', binary_in))
        f.write(namelist('bse_in', bse_settings))
        f.write(namelist('stellar_engine_in', stellar_engine))
    print(f"Wrote {os.path.join(out_dir, 'binary.in')} for bin_num={r.get('bin_num', bin_num)}")

    tracks_file = os.path.join(out_dir, 'metisse_tracks.bin')
    if stellar_engine['using_metisse']:
        recorder, m_min, m_max = record_metisse_tracks(stellar_engine, binary_in['z'])
        write_metisse_tracks(recorder, tracks_file)
        print(f"Wrote {tracks_file} (hydrogen tracks span {m_min}-{m_max} Msun"
              f"{', with helium tracks' if True in recorder.tracks else ', no helium tracks'})")
    elif os.path.exists(tracks_file):
        # avoid confusion with tracks left over from a previous METISSE binary
        os.remove(tracks_file)


def main():
    parser = argparse.ArgumentParser(description='Convert a row from an initC file to a binary.in file')
    parser.add_argument('-f', '--initC_file', type=str, required=True,
                        help='Path to the initC file (hdf5 with key "initC", or csv)')
    parser.add_argument('-b', '--bin_num', type=int, default=None,
                        help='The binary number to convert (optional if there is only one binary)')
    parser.add_argument('-o', '--out_dir', type=str, default=os.path.dirname(os.path.abspath(__file__)),
                        help='Where to write binary.in (default: the debug directory)')
    args = parser.parse_args()

    convert_initC_row_to_binary_in(args.initC_file, args.bin_num, args.out_dir)


if __name__ == '__main__':
    main()
