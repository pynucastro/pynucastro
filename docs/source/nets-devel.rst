****************************
Working with Network Objects
****************************

Rate lists
==========

At its core, a :py:obj:`RateCollection
<pynucastro.networks.rate_collection.RateCollection>` holds a ``list``
of ``Rate`` objects, and uses these to determine everything about the
network.

These are the most important:

* ``RateCollection.rates`` : this is the list of all of the rates in the network
  that *explicitly* link together nuclei in the network.

  This means that rates that are only used as part of rate approximations
  in :py:obj:`ApproximateRate <pynucastro.rates.approximate_rates.ApproximateRate>`
  or as the original / underlying rate in :py:obj:`ModifiedRate <pynucastro.rates.modified_rate.ModifiedRate>` or :py:obj:`BranchedRate <pynucastro.rates.branched_rate.BranchedRate>`
  are not included.

  .. note::

     The method :py:meth:`RateCollection.get_rates <pynucastro.networks.rate_collection.RateCollection.get_rates>` returns this list.

* ``RateCollection.all_rates`` : this is the list of every rate that will be evaluated
  when integrating the network, including those that are *hidden* (i.e., part of ``ApproximateRate``, ``ModifiedRate``, or ``BranchedRate``).

  The rates that are *hidden*, and therefore do not appear in ``RateCollection.rates``
  are identified by ``RateCollection._classify_hidden_rate``, and they is given the 
  ``Rate.removed = True`` attribute.

  It is always that case that:

  .. code:: python

     len(all_rates) >= len(rates)

  .. note::

     The list of rates that are hidden, and therefore not contained in
     ``RateCollection.rates`` can be obtained via
     :py:meth:`RateCollection.get_hidden_rates
     <pynucastro.networks.rate_collection.RateCollection.get_hidden_rates>`.

* type-specific lists (``reaclib_rates``, ``starlib_rates``, ``tabular_rates``, ...):
  these contain just the rates of a particular type or ``Rate`` subclass.  They are
  sorted into these lists by ``RateCollection._build_collection()``.


Updating the network meta-data
==============================

Anytime rates are added, nuclei are added, or other fundamental changes are
done to the network, we need to rerun ``RateCollection._build_collection`` to
update all of the meta-data.

Most functions already handle this themselves.


Ordering in which we fill rates
===============================

In the righthand side function that is output by the different backends, we evaluate
the rates in the following order:

* :py:obj:`ReacLibRate <pynucastro.rates.reaclib_rate.ReacLibRate>`
* :py:obj:`TabularWeakRate <pynucastro.rates.tabular_rate.TabularWeakRate>`
* :py:obj:`TemperatureTabularRate <pynucastro.rates.temperature_tabular_rate.TemperatureTabularRate>`
  and :py:obj:`StarLibRate <pynucastro.rates.starlib_rate.StarLibRate>`
* custom rates (in python only)
* :py:obj:`ModifiedRate <pynucastro.rates.modified_rate.ModifiedRate>`
* :py:obj:`BranchedRate <pynucastro.rates.branched_rate.BranchedRate>`
* :py:obj:`DerivedRate <pynucastro.rates.derived_rate.DerivedRate>`
* :py:obj:`ApproximateRate <pynucastro.rates.approximate_rates.ApproximateRate>`

.. important::

   From this ordering, it is clear that only an ``ApproximateRate`` can depend on
   a ``DerivedRate``.
