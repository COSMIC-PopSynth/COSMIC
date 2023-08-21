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
    'src/dgcore.f',
    'src/evolv2.f',
    'src/gntage.f',
    'src/instar.f',
    'src/kick.f',
    'src/mix.f',
    'src/mrenv.f',
    'src/ran3.f',
    'src/rl.f',
    'src/concatkstars.f',
    'src/comprad.f',
    'src/bpp_array.f',
    'src/checkstate.f',
    'src/deltat.f',
    'src/mlwind.f',
    'src/hrdiag.f',
    'src/star.f',
    'src/zcnsts.f',
    'src/assign_remnant.f',
    'src/hrdiag_remnant.f',
    'src/SSE/SSE_deltat.f',
    'src/SSE/SSE_mlwind.f',
    'src/SSE/SSE_hrdiag.f',
    'src/SSE/SSE_star.f',
    'src/SSE/SSE_zcnsts.f',
    'src/SSE/SSE_zfuncs.f',
    'src/METISSE/src/METISSE_deltat.f90',
    'src/METISSE/src/METISSE_mlwind.f90',
    'src/METISSE/src/METISSE_hrdiag.f90',
    'src/METISSE/src/METISSE_star.f90',
    'src/METISSE/src/METISSE_zcnsts.f90',
    'src/METISSE/src/track_support.f90',
    'src/METISSE/src/z_support.f90',
    'src/METISSE/src/sse_support.f90', 
    'src/METISSE/src/remnant_support.f90',
    'src/METISSE/src/interp_support.f90',
    'src/METISSE/src/alloc_track.f90',
    'src/METISSE/src/dealloc_track.f90',
    'src/METISSE/src/assign_commons_main.f90',
    'src/METISSE/src/assign_commons_COSMIC.f90',
    'src/METISSE/src/comenv_lambda.f90',
    'src/METISSE/src/initialize_front_end.f90']

py3.extension_module('_evolvebin',
  f2py_source,
  lib_source,
  incdir_f2py / 'fortranobject.c',
  include_directories: inc_np,
  install : true,
  subdir : 'cosmic'
)

#fortran_args : ['-Wall']
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
