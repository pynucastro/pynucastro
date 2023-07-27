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
  - For AMReX-Astro Microphysics networks, `examples/triple-alpha/triple-alpha-cxx.py <https://github.com/pynucastro/pynucastro/blob/main/examples/triple-alpha/triple-alpha-cxx.py>`_

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

The `examples/triple-alpha/triple-alpha-cxx.py
<https://github.com/pynucastro/pynucastro/blob/main/examples/triple-alpha/triple-alpha-cxx.py>`_
example builds the right hand side, Jacobian, and helper C++ modules
to copy into the ``networks/`` subdirectory of the AMReX-Astro
Microphysics repository.

No additional customization is required after running the steps for
all networks above.

