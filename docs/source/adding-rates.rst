********************************
Adding Rates from the Literature
********************************

New rate measurements that appear in papers can take a while (years)
before they are picked up by one of the major libraries (ReacLib or
StarLib).

If these rates are given in ReacLib fit format or as a tabulation
(temperature, rate), then it is very easy to add these to pynucastro
and use them in your networks.

Rates in ReacLib format
=======================

If a paper provides a rate in ReacLib format, it is straightforward to
add the rate following the model in :py:mod:`alternate_rates
<pynucastro.rates.alternate_rates>`, but subclassing
:py:obj:`ReacLibRate <pynucastro.rates.reaclib_rate.ReacLibRate>`.
Then you just need to pass the reactants, products, and
:py:obj:`SingleSet <pynucastro.rates.reaclib_rate.SingleSet>` objects
into the constructor.  An example is, for a hypothetical
process $A(B,\gamma)C$ would be:

.. code:: python

   class NewRate(ReacLibRate):

      def __init__(self):

          reactants = ["A", "B"]
          products = ["C"]

          labelprops = "newrate"

          sets = []
          sets.append(SingleSet([a0, a1, a2, a3, a4, a5, a6],
                                labelprops=labelprops))

          super().__init__(reactants=reactants,
                           products=products,
                           sets=sets,
                           label="newrate")

Where ``sets`` can contain as many ``SingleSet`` objects as needed to
describe the rate.

Tabulated Rates
===============

For rates that are given in terms of $(T, \lambda)$ pairs, we can
subclass :py:obj:`TemperatureTabularRate
<pynucastro.rates.temperature_tabular_rate.TemperatureTabularRate>`.

.. important::

   The data must be given as natural logs, e.g., $(\ln(T_9),
   \ln(\lambda))$.

For example, for a hypothetical process $A(B,\gamma)C$ we could do:

.. code:: python

   class NewRate(TemperatureTabularRate):

      def __init__(self):

          reactants = ["A", "B"]
          products = ["C"]

          labelprops = "newrate"

          T9 = np.array([1, 2, 3, 4])
          rate = np.array([1.e-20, 1.e-15, 1.e-10, 1.e-5])

          super().__init__(np.log(T9), np.log(rate),
                           reactants=reactants,
                           products=products,
                           label="newrate")

where the rate data is provided in the NumPy arrays ``T9`` and ``rate``.
