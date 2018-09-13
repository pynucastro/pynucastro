.. pynucastro documentation master file, created by
   sphinx-quickstart on Wed Dec 27 11:54:03 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

pynucastro
==============

`http://github.com/pynucastro/pynucastro <http://github.com/pynucastro/pynucastro>`_

.. note::

   pynucastro does not yet calculate nuclear partition functions to
   correct reverse rates in the Reaclib library. Work is ongoing to
   implement this functionality. We recommend you consider what
   problem you wish to study using pynucastro to determine whether
   reverse rates and partition function corrections are significant at
   the temperatures of interest.

pynucastro is a set of python interfaces to nuclear reaction rate
databases (including the JINA reaclib nuclear reaction rate database).
It is meant for both interactive exploration of rates (through Jupyter
notebooks) and to create reaction networks for use in simulation
codes.


.. toctree::
   :maxdepth: 2
   :caption: Contents

   intro
   pynucastro-examples.ipynb
   library-examples-nuclei.ipynb
   library-examples-filtering.ipynb
   networks
   API <modules>
   sources

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
