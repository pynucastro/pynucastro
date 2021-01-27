Creating Networks
=================

pynucastro can output the righthand side functions for ODE integration
in Python and Fortran.

For Fortran code generation, pynucastro offers two options.

- Standalone: Generate right hand side routines and an integration
  driver program using the included VODE integrator in
  ``util/VODE/``. These networks implement only species and energy
  integration at constant temperature, as pynucastro does not include
  an equation of state.

- Microphysics: Generate network files to copy directly into the
  `StarKiller Microphysics repository <https://github.com/StarKiller-astro/Microphysics/>`_ for use with astrophysical
  simulation codes. These networks implement species, temperature, and
  energy integration.

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
  - For standalone Fortran networks, `examples/urca-23_f90/urca.py <https://github.com/pynucastro/pynucastro/blob/main/examples/urca-23_f90/urca.py>`_
  - For StarKiller Microphysics networks, `examples/urca-23_starkiller/urca.py <https://github.com/pynucastro/pynucastro/blob/main/examples/urca-23_starkiller/urca.py>`_

* Run your python script (pynucastro must be in your PYTHONPATH).

  .. code-block:: bash

     $ python mynet.py

Python network
--------------

First follow the steps above for all networks. To integrate a Python
network using the SciPy integration routines, customize
`examples/CNO/burn.py <https://github.com/pynucastro/pynucastro/blob/main/examples/CNO/burn.py>`_ to initialize and run your network using the
right hand side module you generated above.

Standalone Fortran network
--------------------------

The `examples/urca-23_f90
<https://github.com/pynucastro/pynucastro/tree/main/examples/urca-23_f90>`_
example builds a Fortran 90 network together with a GNU Makefile and
integration driver program using the included VODE package for ODE
integration.

* In the steps for all networks above, pynucastro will create several
  Fortran 90 files as well as a GNU Makefile set up to compile the
  integrator using gfortran. So next do:

  .. code-block:: bash

     $ make

* Edit the generated ``inputs_integration`` file to specify your initial
  conditions as well as integration tolerances, stop time, and output
  sampling interval.

* Run the integrator as:

  .. code-block:: bash

     $ ./main.Linux.gfortran.exe inputs_integration

* The integration results will be stored in a text file named by
  default ``integration_history.dat``

Notes on the build system:

* If you would like to compile with debugging symbols, do:

  .. code-block:: bash

	 $ make NDEBUG=

* To remove the executable and object files produced during compilation, do:

  .. code-block:: bash

     $ make realclean

StarKiller Microphysics network
-------------------------------

The `examples/urca-23_starkiller <https://github.com/pynucastro/pynucastro/tree/main/examples/urca-23_starkiller>`_ example builds the right hand side, Jacobian,
and helper Fortran modules to copy into the ``networks/`` subdirectory
of the StarKiller Microphysics repository.

No additional customization is required after running the steps for
all networks above.

Tabular Rates
-------------

Tabular rates for reactions of the form :math:`\rm{A \rightarrow B}`
are supported by the standalone Fortran and StarKiller Microphysics
network outputs.

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
