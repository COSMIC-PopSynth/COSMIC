project('cosmic', 
    'c',
    'fortran',
    version : '3.4.10',
    default_options: ['warning_level=0', 'optimization=3'],
)

# Enable fortran and check arguments
# add_languages('fortran', native: false)
ff = meson.get_compiler('fortran')
f_args = ff.get_supported_arguments('-fPIC')
add_project_arguments(f_args, language: 'fortran')

py3 = import('python').find_installation()
message(py3.path())


numpy_include_dir = run_command(py3, ['-c', 'import numpy; print(numpy.get_include())'], check: true).stdout().strip()
f2py_include_dir = run_command(py3, ['-c', 'import numpy.f2py; print(numpy.f2py.get_include())'], check: true).stdout().strip()
inc_np = include_directories(numpy_include_dir, f2py_include_dir)

f2py_source = custom_target(
    'evolvebin-target',
    input : ['src/cosmic/src/evolv2.f', 'src/cosmic/src/comprad.f'],
    output : ['_evolvebinmodule.c', '_evolvebin-f2pywrappers.f'],
    command : [py3, '-m', 'numpy.f2py', '@INPUT@', '-m', '_evolvebin', '--lower', '@OUTDIR']
)

lib_source = [
    'src/cosmic/src/hrdiag_remnant.f',
    'src/cosmic/src/assign_remnant.f',
    'src/cosmic/src/benchmarkevolv2.f',
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
    'src/cosmic/src/deltat.f',
    'src/cosmic/src/mlwind.f',
    'src/cosmic/src/hrdiag.f',
    'src/cosmic/src/star.f',
    'src/cosmic/src/zcnsts.f',
    'src/cosmic/src/zfuncs.f',]

# Detect operating system and set appropriate linker flags
host_system = host_machine.system()

if host_system == 'darwin'
  ldflags = ['-Wl,-no_compact_unwind']
else
  ldflags = []  # No special flags for other systems
endif

py3.extension_module('_evolvebin',
  f2py_source,
  lib_source,
  f2py_include_dir / 'fortranobject.c',
  include_directories: inc_np,
  link_args: ldflags,
  install : true,
  install_dir : py3.get_install_dir() / 'cosmic'
)


module_dirs = ['src/cosmic', 'src/cosmic/bse_utils', 
               'src/cosmic/sample']

# Install modules
foreach mod_dir: module_dirs
    install_subdir(mod_dir,
                  install_dir: py3.get_install_dir())
endforeach

message('bin', get_option('bindir'))
python_script = 'bin/cosmic-pop'
install_data(python_script, install_dir: get_option('bindir'))
