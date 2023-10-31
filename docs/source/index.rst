.. pynucastro documentation main file, created by
   sphinx-quickstart on Wed Dec 27 11:54:03 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

pynucastro
==========

`https://github.com/pynucastro/pynucastro <https://github.com/pynucastro/pynucastro>`_

pynucastro is a set of python interfaces to nuclear reaction rate
databases (including the JINA reaclib nuclear reaction rate database).
It is meant for both interactive exploration of rates (through Jupyter
notebooks) and to create reaction networks for use in simulation
codes.


.. toctree::
   :maxdepth: 1
   :caption: Introduction

   intro
   data_sources
   install

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
   electron-capture-example.ipynb
   electron-captures.ipynb

.. toctree::
   :maxdepth: 1
   :caption: Advanced Usage

   screening-examples
   modify-example.ipynb
   approx-rates-examples.ipynb
   custom-rates.ipynb
   unimportant-rates.ipynb
   partition-function

.. toctree::
   :maxdepth: 1
   :caption: Nuclear Statistical Equilibrium

   NSE-example.ipynb

.. toctree::
   :maxdepth: 1
   :caption: Writing and Integrating Networks

   networks
   cxx-networks
   julia

.. toctree::
   :maxdepth: 1
   :caption: Reference

   citing
   API <modules>
   sources

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
