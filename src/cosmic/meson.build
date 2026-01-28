python_sources = [
  '__init__.py',
  '_commit_hash.py',
  '_version.py',
  'checkstate.py',
  'evolve.py',
  'filter.py',
  'get_commit_hash.py',
  'Match.py',
  'plotting.py',
  'utils.py',
  'output.py',
]

py3.install_sources(
  python_sources,
  subdir: 'cosmic'
)

install_data(
  'data/cosmic-settings.json',
  install_dir: py3.get_install_dir() / 'cosmic' / 'data',
)

subdir('bse_utils')
subdir('sample')
subdir('tests')
