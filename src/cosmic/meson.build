numpy = dependency('numpy', required: true)

#Include directories
if not numpy.found()
  error('numpy is required but was not found')
endif

incdir_numpy = run_command(py3,
  ['-c', 'import os; import numpy; print(numpy.get_include())'],
  check: true
).stdout().strip()

incdir_f2py = run_command(py3,
  ['-c', 'import numpy.f2py; print(numpy.f2py.get_include())'],
  check: true
).stdout().strip()
inc_np = include_directories(incdir_numpy, incdir_f2py)

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
  incdir_f2py / 'fortranobject.c',
  include_directories: inc_np,
  link_args: ldflags,
  install : true,
  install_dir: py3.get_install_dir() / 'cosmic'
)

