#Include directories
incdir_numpy = run_command(py3,
  ['-c', 'import os; os.chdir(".."); import numpy; print(numpy.get_include())'],
  check : true
).stdout().strip()

incdir_f2py = run_command(py3,
    ['-c', 'import numpy.f2py; print(numpy.f2py.get_include())'],
    check : true
).stdout().strip()

inc_np = include_directories(incdir_numpy, incdir_f2py)

f2py_source = custom_target(
    'evolvebin-target',
    input : ['src/evolv2.f'],
    output : ['_evolvebinmodule.c', '_evolvebin-f2pywrappers.f'],
    command : [py3, '-m', 'numpy.f2py', '@INPUT@', '-m', '_evolvebin', '--lower', '--build-dir', 'src/cosmic', '@OUTDIR']
)

lib_source = [
    'src/benchmarkevolv2.f',
    'src/corerd.f',
    'src/comenv.f',
    'src/deltat.f',
    'src/dgcore.f',
    'src/evolv2.f',
    'src/gntage.f',
    'src/hrdiag.f',
    'src/instar.f',
    'src/kick.f',
    'src/mix.f',
    'src/mlwind.f',
    'src/mrenv.f',
    'src/ran3.f',
    'src/rl.f',
    'src/star.f',
    'src/zcnsts.f',
    'src/zfuncs.f',
    'src/concatkstars.f',
    'src/comprad.f',
    'src/bpp_array.f',
    'src/checkstate.f',
    'METISSE/src/track_support.f90',
    'METISSE/src/z_support.f90',
    'METISSE/src/zcnsts_METISSE.f90'
]
py3.extension_module('_evolvebin',
  f2py_source,
  lib_source,
  incdir_f2py / 'fortranobject.c',
  include_directories: inc_np,
  install : true,
  subdir : 'cosmic',
  fortran_args : ['-Wall']
)

# Now the source files
py3.install_sources(
    [
        '__init__.py',
        '_version.py',
        'checkstate.py',
        'evolve.py',
        'filter.py',
        'Match.py',
        'plotting.py',
        'utils.py',
        'bse_utils/zcnsts.py',
        'bse_utils/__init__.py',
        'bse_utils/zdata.py'
    ],
    pure: false,
    subdir: 'cosmic'
)
