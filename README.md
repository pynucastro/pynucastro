# pynucastro

[![Build Status](https://travis-ci.org/pynucastro/pynucastro.svg?branch=master)](https://travis-ci.org/pynucastro/pynucastro) [![JOSS Status](http://joss.theoj.org/papers/f753b6f21f460ae6a301c21c95dfa001/status.svg)](http://joss.theoj.org/papers/f753b6f21f460ae6a301c21c95dfa001)

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

  * [pynucastro-examples.ipynb](https://github.com/pynucastro/pynucastro/blob/master/examples/pynucastro-examples.ipynb)

  * [interactive-example.ipynb](https://github.com/pynucastro/pynucastro/blob/master/examples/interactive-example.ipynb)


# setup

To set this up, you should include the root pynucastro path in your
`PYTHONPATH` environment variable.


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
py.test --nbval examples/interactive-example.ipynb examples/pynucastro-examples.ipynb examples/pp-CNO-example.ipynb
```

If your OS has both Python 2 and Python 3 installed you may need to
invoke `pytest` as `py.test-3` when running the unit tests.
