Common Networks
===============

pynucastro provides the :py:mod:`common_networks <pynucastro.networks.common_networks>` module
to automatically create some commonly-used networks.

The current networks are:

* ``aprox13`` : the 13-isotope $\alpha$-chain network.  This can be created via:

  .. code:: python

     net = pyna.common_networks.aprox13()

  See :doc:`aprox13` for details on the approximations used in this network.

* MESA's ``basic.net`` : the 8-isotope H/He burning network that MESA uses for
  basic stellar evolution.  This can be created via:

  .. code:: python

     net = pyna.common_networks.mesa_basic()

  See :doc:`mesa-basic` for details on the approximations used in this network.

* ``cno`` : a general network appropriate to CNO/hot-CNO burning with rp-breakout.
  This is based on the MESA ``cno_extras`` network, but it has a few more reverse
  rates.

  .. code:: python

     net = pyna.common_networks.cno()

  See :doc:`cno-network` for details on the approximations used in this network.

