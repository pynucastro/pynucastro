# pynucastro

[![Build Status](https://travis-ci.org/pynucastro/pynucastro.svg?branch=main)](https://travis-ci.org/pynucastro/pynucastro)
[![DOI](http://joss.theoj.org/papers/10.21105/joss.00588/status.svg)](https://doi.org/10.21105/joss.00588)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1202434.svg)](https://doi.org/10.5281/zenodo.1202434)
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/pynucastro/pynucastro/main?filepath=examples%2Fpynucastro-examples.ipynb)

Python interfaces to the nuclear reaction rate databases, including
the JINA Reaclib nuclear reactions database.  This
project has 2 goals:

  * allow for an exploration of rates and collection of rates (networks)
    in Jupyter notebooks

  * allow for the easy creation of the righthand side routines for
    reaction network integration (the ODEs) 

To do this, pynucastro provides a parser for the Reaclib format to
produce either a callable python function for a rate or output the
python code for a function that can be incorporated into a rate
module.

pynucastro does not yet support nuclear partition functions for
Reaclib reverse rates, as the implementation is currently under
development. We recommend you consider what problem you wish to study
using pynucastro to determine whether reverse rates and partition
function corrections are significant at the temperatures of interest.

# documentation

Documentation for pynucastro is available here:

http://pynucastro.github.io/pynucastro/


# example
```
import pynucastro

r = pynucastro.Rate("c13-pg-n14-nacr")

# evaluate this rate at T = 1.e9 K
print(r.eval(1.0e9))

```

Also see the Jupyter notebooks:

  * [pynucastro-examples.ipynb](https://github.com/pynucastro/pynucastro/blob/main/examples/pynucastro-examples.ipynb)

  * [interactive-example.ipynb](https://github.com/pynucastro/pynucastro/blob/main/examples/interactive-example.ipynb)


# install

To install the package, you can run:
```
python setup.py install
```
for a systemwide install, or
```
python setup.py install --user
```
for a single-user install.  This will put the pynucastro modules and library in
the default location python searches for packages.


# requirements

This package requires Python 3 (release 3.4 or later) and the following Python packages:

* `numpy`

* `sympy`

* `scipy`

* `matplotlib`

* `networkx`

* `ipywidgets`

To build the documentation or run the unit tests, `sphinx` and
`pytest` are additionally required along with some supporting
packages. See the included `requirements.txt` file for a list of these
packages and versions. To install the packages from the requirements
file, do:
```
pip install -r requirements.txt
```

# unit tests

We use py.test to do unit tests.  In `pynucastro/`, do:
```
py.test -v .
```

to see coverage, do:
```
py.test --cov=pynucastro .
```

to test the notebooks, do:
```
py.test --nbval examples
```

If your OS has both Python 2 and Python 3 installed you may need to
invoke `pytest` as `py.test-3` when running the unit tests.


# Core Developers

People who make a number of substantive contributions (new features,
bug fixes, etc.) will be named "core developers" of pynucastro.

Core developers will be recognized in the following ways:

  * invited to the group's slack team

  * listed in the User's Guide and website as a core developer

  * listed in the author list on the Zenodo DOI for the project
    (as given in the .zenodo.json file)

  * invited to co-author general code papers / proceedings describing
    pynucastro.  (Note: science papers that use pynucastro will always
    be left to the science paper lead author to determine authorship).

If a core developer is inactive for 3 years, we may reassess their
status as a core developer.

