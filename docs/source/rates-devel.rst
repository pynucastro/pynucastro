**************************************
Notes on Working with ``Rate`` Objects
**************************************

Attributes
==========

There are several different attributes that refer to the name / properties of the rate in some way.  These are:

* ``Rate.id``

* ``Rate.rid``

* ``Rate.fname``


Copying
=======

A ``Rate`` is not immutable.  This can be an issue when we we do:

.. code:: python

   rl = pyna.ReacLibLibrary()
   c12ag = rl.get_rate_by_name("c12(a,g)o16")

Here, the rate ``c12ag`` is actually a reference to the version stored
internally in the ``Library`` that ``ReacLibLibrary`` creates.

This means if we modify the rate, e.g, as:

.. code:: python

   c12ag.removed = True

then our library ``rl`` is also updated.  If we reuse this library,
this can have unintended consequences.

As a result, when working with rates that you think may be modified,
it is best to work with a copy.  The recommended copy is to use
the python ``copy`` module's shallow copy:

.. code:: python

   c12ag_copy = copy.copy(c12ag)

This will invoke a custom ``__copy__`` method in the rate that does a
shallow copy for most rate attributes, but explicitly recreates the
``reactants``, ``products`` and ``stoichiometry`` data.

.. tip::

   When we create a ``RateCollection`` or other network, we explicitly
   copy the rates into the network.
