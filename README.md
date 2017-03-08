# pyreaclib

Python interfaces to the JINA reaclib nuclear reactions database.  This
project has 2 goals:

  * allow for an exploration of rates and collection of rates (networks)
    in Jupyter notebooks

  * allow for the easy creation of the righthand side routines for
    reaction network integration (the ODEs) 

To do this, pyreaclib providea a parser for the reaclib format to
produce either a callable python function for a rate or output the
python code for a function that can be incorporated into a rate
module.



# example

```
import pyreaclib

r = pyreaclib.Rate("c13-pg-n14-nacr")

# evaluate this rate at T = 1.e9 K
print(r.eval(1.0e9))

```

Also see the Jupyter notebooks:

  * [pyreaclib-examples.ipynb](https://github.com/zingale/pyreaclib/blob/master/pyreaclib-examples.ipynb)

  * [interactive-examples.ipynb](https://github.com/zingale/pyreaclib/blob/master/interactive-examples.ipynb)


# setup

To set this up, you should include the root pyreaclib path in your
`PYTHONPATH` environment variable.


# requirements

This package requires

* `numpy`

* `scipy`: only needed for the example integrator

* `periodictable`: needed to get Z from element symbol, necessary for
   sorting

* `matplotlib`: needed for plotting the rate

* `networkx`: for visualizing the network


# TODO

* output Cython code

* return a polynomial fit to a rate in a given interval

* create a Ydot() class that can both output dYdt[nucleus] and
  evaluate dYdt.
  

