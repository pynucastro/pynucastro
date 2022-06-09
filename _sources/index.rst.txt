.. pynucastro documentation main file, created by
   sphinx-quickstart on Wed Dec 27 11:54:03 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

pynucastro
==========

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
   :maxdepth: 1
   :caption: Introduction

   intro

.. toctree::
   :maxdepth: 1
   :caption: Exploring Networks in Python

   pynucastro-examples.ipynb
   plot-types.ipynb
   pynucastro-integration.ipynb

.. toctree::
   :maxdepth: 1
   :caption: Selecting Rates

   library-examples.ipynb
   validate-example.ipynb
   electron-captures.ipynb

.. toctree::
   :maxdepth: 1
   :caption: Advanced Usage

   modify-example.ipynb
   approx-rates-examples.ipynb
   partition-function

.. toctree::
   :maxdepth: 1
   :caption: Writing Networks

   networks
   cxx-networks

.. toctree::
   :maxdepth: 1
   :caption: Reference

   API <modules>
   sources

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
