Installing ``diffeqpy``
=======================

.. important::

   In order to use ``diffeqpy``, julia must be installed and in the path, along with ``DifferentialEquations.jl`` and ``PyCall.jl``

From the julia prompt, we can install the needed packages:

.. code:: Julia

   using Pkg
   Pkg.add("DifferentialEquations")
   Pkg.add("PyCall")`

We can then install ``diffeqpy`` via pip:

.. prompt:: bash

   pip install diffeqpy

It is also suggested that you have Numba installed for better performance.  Numba is a dependency
of pynucastro so it is likely already installed.

