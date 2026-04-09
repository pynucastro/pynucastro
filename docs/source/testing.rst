******************************
Testing and Comparing Backends
******************************

If we create a network with a set of rates, we should get the same
"ydots" ($dY_k/dt$) values regardless of the backend (to roundoff, and maybe
interpolation, error)

To facilitate comparing the backends, we have a ``NetworkCompaare``
helper class that takes a library and evaluates the ydots for different backends.
Currently it can test:

* inline python networks: this is what we do with a ``RateCollection``
  using the built-in ``eval`` methods to evaluate each of the rates.

* ``PythonNetwork`` modules: this means writing out the network as a
  ``.py`` file, importing it, and then using the functions in the
  module to evaluate the rates.

* ``SimpleCxxNetwork`` : we build a simple test driver and parse the output
  to get the ydot values.

This is used in the following unit tests in ``pynucastro/networks/tests/``:

* ``test_compare_cxx_and_python.py`` : this compares just strong reaction rates

* ``test_compare_temperature_tabular.py`` : this compares the evaluation of
  a ``TemperatureTabularRate`` between the inline and module python nets.

As more features are ported to ``SimpleCxxNetwork``, we should extend the
testing to compare with python rates.

