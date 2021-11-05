import re
from setuptools import setup, find_packages

VERSIONFILE = "pynucastro/_version.py"
VERSION_RE = r"^__version__ = ['\"]([^'\"]*)['\"]"

def get_version():
    try:
        with open(VERSIONFILE) as vf:
            line = vf.readline()

        version_group = re.search(VERSION_RE, line, re.M)
        return version_group.group(1)
    except FileNotFoundError:
        # this can be triggered by the docs build
        return ""

setup(name='pynucastro',
      version=get_version(),
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
      zip_safe=False)
