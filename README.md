# pyreaclib

Python interfaces to the JINA reaclib nuclear reactions database.  The
main goal is to provide a parser for the reaclib format to produce
either a callable python function for a rate or output the python code
for a function that can be incorporated into a rate module.

For performance reasons, the common temperature factors needed by a
rate are stored in a `Tfactors` class---a `Tfactors` object is passed
into each rate function to get the rate.


# example

```
import reaclib

r = reaclib.Rate("c13-pg-n14-nacr")

# output python code for this rate
print r.function_string()

# evaluate this rate at T = 1.e9 K
print r.eval(1.0e9)

```

Also see the Jupyter notebook:
[pyreaclib-examples.ipynb](https://github.com/zingale/pyreaclib/blob/master/pyreaclib-examples.ipynb)


# requirements

This package requires

* `numpy`

* `scipy`: only needed for the example integrator

* `periodictable`: needed to get Z from element symbol, necessary for
   sorting

* `matplotlib`: needed for plotting the rate

* `networkx`: for visualizing the network


# TODO

* output Fortran functions as an option

  - integrate OpenACC support

* output Cython code

* return a polynomial fit to a rate in a given interval

* create a Ydot() class that can both output dYdt[nucleus] and
  evaluate dYdt.

* draw a diagram of the nuclei in a RateCollection along with arrows
  indicating the flows

  - use NetworkX for this -- nuclei are nodes, positioned along
    Z-N plane and colored according to dYdt


