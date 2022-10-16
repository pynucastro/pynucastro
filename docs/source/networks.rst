Creating Networks
=================

pynucastro can output the righthand side functions for ODE integration
in Python or in a form for the AMReX-Astrophysics Microphysics package in C++


All networks
------------

The first few steps to setting up any of the available network options
are as follows. See the following section for details about the
particular kind of network you are interested in creating.

* Create a working directory for your network (this does not have to
  be in the ``pynucastro/`` directory tree).

* Obtain Reaclib rate files (in Reaclib 1 or 2 format) for your problem and
  either place them in your working directory or ``pynucastro/library``.

  You can find the rates that are already provided in the library by
  doing:

  .. code-block:: python

     pynucastro.list_known_rates()


* Write a short python script to generate your network,
  e.g. ``mynet.py``, based on the following example scripts:

  - For Python networks, `examples/CNO/cno.py <https://github.com/pynucastro/pynucastro/blob/main/examples/CNO/cno.py>`_
  - For AMReX-Astro Microphysics networks, `examples/urca-23_starkiller/urca.py <https://github.com/pynucastro/pynucastro/blob/main/examples/urca-23_starkiller/urca.py>`_

* Run your python script

  .. code-block:: bash

     $ python mynet.py

Python network
--------------

First follow the steps above for all networks. To integrate a Python
network using the SciPy integration routines, customize
`examples/CNO/burn.py <https://github.com/pynucastro/pynucastro/blob/main/examples/CNO/burn.py>`_ to initialize and run your network using the
right hand side module you generated above.


AMReX-Astro Microphysics network
--------------------------------

The `examples/urca-23_starkiller <https://github.com/pynucastro/pynucastro/tree/main/examples/urca-23_starkiller>`_ example builds the right hand side, Jacobian,
and helper C++ modules to copy into the ``networks/`` subdirectory
of the AMReX-Astro Microphysics repository.

No additional customization is required after running the steps for
all networks above.

Tabular Rates
-------------

Tabular rates for reactions of the form :math:`\rm{A \rightarrow B}`
are supported by the AMReX-Astro Microphysics network outputs.

If you would like to include tabular rates, for now they must be in
the form of, e.g. `23Na-23Ne_electroncapture.dat <https://github.com/pynucastro/pynucastro/blob/main/pynucastro/library/tabular/23Na-23Ne_electroncapture.dat>`_ in
``pynucastro/library/tabular/``, indexed by the product of density and
electron fraction :math:`\rm{\rho Y_e}` and temperature
:math:`\rm{T}`, with the same number and order of variables.

To generate a network with a tabular rate, prepare a rate file to
describe how to read the table as below and then list it as you would
a Reaclib rate file in your network generation script. For example,
`pynucastro/library/tabular/na23--ne23-toki <https://github.com/pynucastro/pynucastro/blob/main/pynucastro/library/tabular/na23--ne23-toki>`_ demonstrates the following
format:

.. code-block:: none

   t
   [parent nuclide]  [daughter nuclide]
   [rate table file name]
   [number of header lines before the first line of data]
   [number of density*ye values]
   [number of temperature values]
