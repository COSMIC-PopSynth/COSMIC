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

py_mod = import('python')
py3 = py_mod.find_installation('python3')

numpy_include_dir = run_command(py3, ['-c', 'import numpy; print(numpy.get_include())'], check: true).stdout().strip()
f2py_include_dir = run_command(py3, ['-c', 'import numpy.f2py; print(numpy.f2py.get_include())'], check: true).stdout().strip()
message('f2py', f2py_include_dir)
inc_np = include_directories(numpy_include_dir, f2py_include_dir)

f2py_source = custom_target(
    'evolvebin-target',
    input : ['src/evolv2.f', 'src/comprad.f'],
    output : ['_evolvebinmodule.c', '_evolvebin-f2pywrappers.f'],
    command : [py3, '-m', 'numpy.f2py', '@INPUT@', '-m', '_evolvebin', '--lower', '--build-dir', 'src/cosmic']
)

lib_source = [
    'src/hrdiag_remnant.f',
    'src/assign_remnant.f',
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
    'src/deltat.f',
    'src/mlwind.f',
    'src/hrdiag.f',
    'src/star.f',
    'src/zcnsts.f',
    'src/zfuncs.f',]

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
  install : true
)


module_dirs = ['src/cosmic', 'src/cosmic/bse_utils', 
               'src/cosmic/sample']

# Install modules
foreach mod_dir: module_dirs
    install_subdir(mod_dir,
                  install_dir: py3.get_install_dir())
endforeach

python_script = 'bin/cosmic-pop'

install_data(python_script, install_dir: get_option('bindir'))