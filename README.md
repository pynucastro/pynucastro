# pynucastro

[![Docs](https://github.com/pynucastro/pynucastro/workflows/github%20pages/badge.svg)](http://pynucastro.github.io/pynucastro/)
[![pytest-all](https://github.com/pynucastro/pynucastro/actions/workflows/pytest-all.yml/badge.svg?branch=main)](https://github.com/pynucastro/pynucastro/actions/workflows/pytest-all.yml)
[![pylint](https://github.com/pynucastro/pynucastro/actions/workflows/pylint.yml/badge.svg?branch=main)](https://github.com/pynucastro/pynucastro/actions/workflows/pylint.yml)
[![flake8](https://github.com/pynucastro/pynucastro/actions/workflows/flake8.yml/badge.svg?branch=main)](https://github.com/pynucastro/pynucastro/actions/workflows/flake8.yml)
[![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)
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

An extensive demonstration of the capabilities of pynucastro is shown in this notebook:

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/pynucastro/pynucastro/main?filepath=examples%2Fpynucastro-examples.ipynb)

[pynucastro-examples.ipynb](https://github.com/pynucastro/pynucastro/blob/main/examples/pynucastro-examples.ipynb)


We can also interactively explore a reaction network.  Here's an example of hot-CNO with breakout reactions:

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/pynucastro/pynucastro/HEAD?labpath=examples%2Fhot-CNO-breakout-example.ipynb)

by increasing the temperature, you can see the transition from CNO to
hot-CNO (proton captures on C13 become faster than the beta decay) and
then the breakout of hot-CNO when the alpha capture on O14 becomes
faster than its decay.


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

This package requires Python 3 (release 3.6 or later) and the following Python packages:

* `numpy`

* `sympy`

* `scipy`

* `matplotlib`

* `networkx`

* `ipywidgets`

* `numba`

* `setuptools_scm`

To build the documentation or run the unit tests, `sphinx` and
`pytest` are additionally required along with some supporting
packages. See the included `requirements.txt` file for a list of these
packages and versions. To install the packages from the requirements
file, do:
```
pip install -r requirements.txt

```
Is important to stress out that all the dependencies should be 
installed before `pynucastro`, otherwise the installation should be 
repeated.

# unit tests

We use py.test to do unit tests.  In `pynucastro/`, do:
```
pytest -v .
```

to see coverage, do:
```
pytest --cov=pynucastro .
```

to test the notebooks, do:
```
py.test --nbval examples
```


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

