from setuptools import setup, find_packages

setup(name='pynucastro',
      version='1.2.0',
      description='Python Interfaces to the nuclear reaction rate databases',
      url='https://github.com/pynucastro/nucastro',
      author='Mike Zingale and Donald Willcox',
      author_email='michael.zingale@stonybrook.edu',
      license='BSD',
      packages=find_packages(),
      package_data={"pynucastro": ["library/*", "library/tabular/*", "templates/*", "templates/fortran-vode/*", "templates/starkiller-microphysics/*", "nucdata/*"]},
      install_requires=['networkx', 'numpy',
                        'sympy', 'scipy', 'matplotlib',
                        'ipywidgets'],
      zip_safe=False)
