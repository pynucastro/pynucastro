Python network
==============

A :class:`PythonNetwork <pynucastro.networks.python_network.PythonNetwork>` will output a python module with all the data
needed for constructing the righthand side and Jacobian of the ODE
system.  This is designed to solve the system:

.. math::

   \frac{dY_i}{dt} = f(\rho, T, \{Y_1, Y_2, ...\})

where $Y_i$ is the molar fraction of species $i$.

.. note::

   The functions output when the network is written are marked up with
   Numba JIT decorators for performance.

Using the previous example:

.. code:: python

   import pynucastro as pyna

   rl = pyna.ReacLibLibrary()
   lib = rl.linking_nuclei(["he4", "c12", "o16"])

   net = pyna.PythonNetwork(libraries=[lib])
   net.write_network("triple_alpha.py")

after running this script, we could then import the network as:

.. code:: python

   import triple_alpha

access the information needed for an ODE integrator.

For this example, 4 rates will be found (forward and reverse).


The following information is provided:

* a set of integer keys for indexing nuclei.  In this case, ``jhe4``,
  ``jc12``, and ``jo16``, along with the number of nuclei, ``nnuc``

* the nuclear data, including arrays of length ``nnuc``:

  * ``A[]`` : the atomic weights
  * ``Z[]`` : the atomic numbers
  * ``mass[]`` : the mass of each nuclei (in ergs)
  * ``names[]`` : the descriptive name (e.g. ``"He4"``) on the nucleus

* some helper functions:

  * ``to_composition(Y)`` : convert an array of molar fractions to a :func:`Composition
    <pynucastro.networks.rate_collection.Composition>` object

  * ``energy_release(dY)`` : computes the energy release (in erg/g) given the change
    in molar fractions ``dY``.

  * ``ye(Y)`` : computes the electron fraction, $Y_e$, for input molar fractions ``Y``.

* rate functions

  To store the rates as they are evaluated, a class called
  ``RateEval`` that has members for each of the rate.  For example, for
  the triple alpha rate in the example above, the rate will be stored as
  ``RateEval.He4_He4_He4__C12``.

  Each rate then has its own function that updates ``RateEval`` directly.  These
  take the form:

  .. code:: python

     rate_name(rate_eval, tf)

  where ``rate_eval`` is a ``RateEval`` object and ``tf`` is a
  :func:`Tfactors <pynucastro.rates.rate.Tfactors>` object

* RHS and Jacobian

  Finally, to interface with an ODE solver, we have:

  * ``rhs(t, Y, rho, T, screen_func=None)`` : the righthand side function.
    This takes time (``t``), the molar fraction array (``Y``), density (``rho``), temperature (``T``),
    and (optionally) any of the screening functions that should be applied.

    It returns the array of $dY/dt$.

  * ``jacobian(t, Y, rho, T, screen_func=None``) : the Jacobian routine.  The arguments are the same
    as for ``rhs``. 

    This returns the Jacobian array, $J_{i,j} = \partial \dot{Y}_i/\partial Y_j$.

A ``PythonNetwork`` can be integrated using the SciPy `solve_ivp <https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html>`_ methods.

