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
py3 = py_mod.find_installation(pure: false)
py3_dep = py3.dependency()
message(py3.path())
message(py3.get_install_dir())

subdir('./src/cosmic')

module_dirs = ['./src/cosmic', './src/cosmic/bse_utils', 
               './src/cosmic/sample']

# Install modules
foreach mod_dir: module_dirs
    install_subdir(mod_dir,
                  install_dir: py3.get_install_dir())
endforeach

python_script = 'bin/cosmic-pop'

install_data(python_script, install_dir: get_option('bindir'))