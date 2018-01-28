# pynucastro

Python interfaces to the nuclear reaction rate databases, including
the JINA reaclib nuclear reactions database.  This
project has 2 goals:

  * allow for an exploration of rates and collection of rates (networks)
    in Jupyter notebooks

  * allow for the easy creation of the righthand side routines for
    reaction network integration (the ODEs) 

To do this, pynucastro providea a parser for the reaclib format to
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

  * [pynucastro-examples.ipynb](https://github.com/zingale/pynucastro/blob/master/pynucastro-examples.ipynb)

  * [interactive-examples.ipynb](https://github.com/zingale/pynucastro/blob/master/interactive-examples.ipynb)


# setup

To set this up, you should include the root pynucastro path in your
`PYTHONPATH` environment variable.


# requirements

This package requires

* `numpy`

* `scipy`: only needed for the example integrator

* `periodictable`: needed to get Z from element symbol, necessary for
   sorting

* `matplotlib`: needed for plotting the rate

* `networkx`: for visualizing the network


# unit tests

We use py.test to do unit tests.  In `pynucastro/`, do:
```
py.test-3 -v .
```

to see coverage, do:
```
py.test-3 --cov=pynucastro .
```

to test the notebooks, do:
```
py.test-3 --nbval examples
```

# TODO

* output Cython code (rhs)

  * probably want to use a cdef class for the TFactors:
    http://stackoverflow.com/questions/21012348/cython-how-to-make-an-python-object-as-a-property-of-cython-class
	http://cython.readthedocs.io/en/latest/src/tutorial/cdef_classes.html
	http://cython.readthedocs.io/en/latest/src/userguide/early_binding_for_speed.html
	



* return a polynomial fit to a rate in a given interval

* create a Ydot() class that can both output dYdt[nucleus] and
  evaluate dYdt.
  

