project('cosmic', 
    'c',
    'fortran',
    version : '3.5.0',
    default_options: ['warning_level=0', 'optimization=3'],
)

# Enable fortran and check arguments
# add_languages('fortran', native: false)
ff = meson.get_compiler('fortran')
f_args = ff.get_supported_arguments('-fPIC')
add_project_arguments(f_args, language: 'fortran')


py3 = import('python').find_installation(pure: false)

numpy_include_dir = run_command(py3, ['-c', 'import numpy; print(numpy.get_include())'], check: true).stdout().strip()
f2py_include_dir = run_command(py3, ['-c', 'import numpy.f2py; print(numpy.f2py.get_include())'], check: true).stdout().strip()
inc_np = include_directories(numpy_include_dir, f2py_include_dir)

get_hash = run_command('python', './src/cosmic/get_commit_hash.py', check: true).stdout().strip()


lib_source = [
    'src/cosmic/src/hrdiag_remnant.f',
    'src/cosmic/src/assign_remnant.f',
    'src/cosmic/src/benchmarkevolv2.f',
    'src/cosmic/src/int64.f',
    'src/cosmic/src/corerd.f',
    'src/cosmic/src/comenv.f',
    'src/cosmic/src/dgcore.f',
    'src/cosmic/src/evolv2.f',
    'src/cosmic/src/gntage.f',
    'src/cosmic/src/instar.f',
    'src/cosmic/src/kick.f',
    'src/cosmic/src/mix.f',
    'src/cosmic/src/mrenv.f',
    'src/cosmic/src/ran3.f',
    'src/cosmic/src/rl.f',
    'src/cosmic/src/concatkstars.f',
    'src/cosmic/src/comprad.f',
    'src/cosmic/src/bpp_array.f',
    'src/cosmic/src/checkstate.f',
    'src/cosmic/src/deltat.f',
    'src/cosmic/src/mlwind.f',
    'src/cosmic/src/hrdiag.f',
    'src/cosmic/src/star.f',
    'src/cosmic/src/zcnsts.f',
    'src/cosmic/src/SSE/SSE_deltat.f',
    'src/cosmic/src/SSE/SSE_mlwind.f',
    'src/cosmic/src/SSE/SSE_hrdiag.f',
    'src/cosmic/src/SSE/SSE_star.f',
    'src/cosmic/src/SSE/SSE_zcnsts.f',
    'src/cosmic/src/SSE/SSE_zfuncs.f',
    'src/cosmic/src/SSE/SSE_gntage.f',
    'src/cosmic/src/METISSE/src/METISSE_gntage.f90',
    'src/cosmic/src/METISSE/src/METISSE_deltat.f90',
    'src/cosmic/src/METISSE/src/METISSE_mlwind.f90',
    'src/cosmic/src/METISSE/src/METISSE_hrdiag.f90',
    'src/cosmic/src/METISSE/src/METISSE_star.f90',
    'src/cosmic/src/METISSE/src/METISSE_zcnsts.f90',
    'src/cosmic/src/METISSE/src/track_support.f90',
    'src/cosmic/src/METISSE/src/z_support.f90',
    'src/cosmic/src/METISSE/src/sse_support.f90', 
    'src/cosmic/src/METISSE/src/remnant_support.f90',
    'src/cosmic/src/METISSE/src/interp_support.f90',
    'src/cosmic/src/METISSE/src/comenv_lambda.f90',
    'src/cosmic/src/METISSE/src/METISSE_miscellaneous.f90',
    'src/cosmic/src/assign_commons_COSMIC.f90']

# Detect operating system and set appropriate linker flags
host_system = host_machine.system()

if host_system == 'darwin'
  ldflags = ['-Wl,-no_compact_unwind']
else
  ldflags = []  # No special flags for other systems
endif


f2py_source = custom_target(
    'evolvebin-target',
    input : ['src/cosmic/src/evolv2.f', 'src/cosmic/src/comprad.f'],
    output : ['_evolvebinmodule.c', '_evolvebin-f2pywrappers.f'],
    command : [py3, '-m', 'numpy.f2py', '@INPUT@', '-m', '_evolvebin', '--lower']
)


evolvebin_module = py3.extension_module('_evolvebin',
  f2py_source,
  lib_source,
  f2py_include_dir / 'fortranobject.c',
  include_directories: inc_np,
  link_args: ldflags,
  install : true,
  subdir : 'cosmic'
)

python_script = 'bin/cosmic-pop'
install_data(python_script, install_dir: get_option('bindir'))

subdir('src/cosmic')
