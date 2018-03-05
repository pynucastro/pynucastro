from setuptools import setup

setup(name='pynucastro',
      version='1.0',
      description='Python Interfaces to the nuclear reaction rate databases',
      url='https://github.com/pynucastro/nucastro',
      author='Mike Zingale and Donald Willcox',
      author_email='michael.zingale@stonybrook.edu',
      license='BSD',
      packages=['pynucastro'],
      install_requires=['networkx', 'numpy',
                        'sympy', 'scipy', 'matplotlib',
                        'ipywidgets']
      zip_safe=False)
