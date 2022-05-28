from setuptools import setup, find_packages

setup(name='pynucastro',
      description='Python Interfaces to nuclear reaction rate databases',
      url='https://github.com/pynucastro/nucastro',
      author='pynucastro development group',
      author_email='michael.zingale@stonybrook.edu',
      license='BSD',
      packages=find_packages(),
      package_data={"pynucastro": ["library/*", "library/tabular/*", "templates/*",
                                   "templates/starkiller-microphysics/*",
                                   "templates/starkiller-cxx-microphysics/*", "nucdata/*",
                                   "nucdata/AtomicMassEvaluation/*", "nucdata/PartitionFunction/*"]},

      install_requires=['networkx', 'numpy', 'sympy',
                        'scipy', 'matplotlib', 'ipywidgets'],
      extras_require={"numba": ["numba"]},
      use_scm_version={"version_scheme": "post-release",
                       "write_to": "pynucastro/_version.py"},
      setup_requires=["setuptools_scm"],
      zip_safe=False)
