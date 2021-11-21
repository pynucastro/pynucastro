from setuptools import setup, find_packages

setup(name='pynucastro',
      version='1.5.0',
      description='Python Interfaces to the nuclear reaction rate databases',
      url='https://github.com/pynucastro/nucastro',
      author='pynucastro development group',
      author_email='michael.zingale@stonybrook.edu',
      license='BSD',
      packages=find_packages(),
      package_data={"pynucastro": ["library/*", "library/tabular/*", "templates/*", "templates/fortran-vode/*", "templates/starkiller-microphysics/*", "nucdata/*"]},
      install_requires=['networkx', 'numpy', 'numba',
                        'sympy', 'scipy', 'matplotlib',
                        'ipywidgets', 'numba'],
      use_scm_version={"version_scheme": "post-release",
                       "write_to": "pynucastro/_version.py"},
      zip_safe=False)
