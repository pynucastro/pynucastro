*****************************
Working with ``Rate`` Objects
*****************************

Descriptive attributes
======================

There are several different attributes that refer to the name /
properties of the rate in some way.  Consider the rate for
${}^{12}\mathrm{C}(\alpha,\gamma){}^{16}\mathrm{O}$:

.. code:: python

   rl = pyna.ReacLibLibrary()
   c12ag = rl.get_rate_by_name("c12(a,g)o16")

The following rate attributes are defined:

* ``Rate.rid`` : this is a string that gives the reactants and
  products, using only ASCII characters.  No information about the
  source of the rate included.

  It is not used for library operations directly, but instead is used
  to build ``rid`` (described next) and included as a human-readable
  comment in the generated rate functions in exported networks.

  For ``c12ag.rid``, we have:

  ::

     'C12 + He4 --> O16'

  .. note::

     Other sources of ${}^{12}\mathrm{C}(\alpha,\gamma){}^{16}\mathrm{O}$ will have the
     same ``rid``.  For instance:

     .. code:: python

        from pynucastro.rates.alternate_rates import DeBoerC12agO16
        c12ag_deboer = DeBoerC12agO16()
      
     We would find ``c12ag_deboer.rid`` to be:

     ::

        'C12 + He4 --> O16'



* ``Rate.id`` : this is intended to be a unique identifier for a rate.
  It is composed of the ``rid`` and the rate source (``Rate.src``).

  Here are several ways ``Rate.id`` is used:

  * To test equality of rates, for example, given a
    ``Rate`` ``r`` and a list of rates ``rate_list`` doing:

    .. code:: python

       if r in rate_list:
           # do something

    will use ``r.id`` for the comparison.

  * For a ``Library``, the main storage for rates (``Library._rates``)
    is a dictionary keyed by ``Rate.id``.

    As a result, ``Rate,id`` must
    be unique in a ``Library`` and a generated network.

  * For a ``RateCollection``, we use lists for ``RateCollection.rates``
    and ``RateCollection.all_rates``.  As noted above, to test membership
    in a list, the ``id`` will be used.

  For ``c12ag.id``, we have:

  ::

     'C12 + He4 --> O16 <reaclib_nac2>'

  This is distinct from other sources of the ${}^{12}\mathrm{C}(\alpha,\gamma){}^{16}\mathrm{O}$ rate,
  so using ``c12ag_deboer`` defined above, we could see ``c12ag_deboer.id`` as:

  ::

     'C12 + He4 --> O16 <deboer_deboer2017>'

* ``Rate.fname`` : this is a programming-language safe name for the
  string (no special characters except for ``_``), and is used for the function names
  and rate indices in exported networks.

  This must be unique in a generated network.

  For ``c12ag``, we have:

  ::

     'He4_C12_to_O16_reaclib'

  For weak rates, the weak rate type is added to the string.  This is important,
  for example, since ReacLib provides two rates for $p + p$, a $\beta^+$ and $e^-$-capture:

  .. code:: python

     pp, pep = rl.get_rate_by_name("p(p,)d")

  If we look at these separately, we see:

  * ``pp.fname`` is ``'p_p_to_d_beta_pos_reaclib'``
  * ``pep.fname`` is ``'p_p_to_d_electron_capture_reaclib'``


Copying
=======

A ``Rate`` is not immutable.  This can be an issue when we do:

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
