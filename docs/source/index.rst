.. pynucastro documentation main file, created by
   sphinx-quickstart on Wed Dec 27 11:54:03 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

pynucastro
==========

.. image:: logo.png

`https://github.com/pynucastro/pynucastro <https://github.com/pynucastro/pynucastro>`_

pynucastro is a set of python interfaces to nuclear reaction rate
databases (including the JINA reaclib nuclear reaction rate database).
It is meant for both interactive exploration of rates (through Jupyter
notebooks) and to create reaction networks for use in simulation
codes.


.. toctree::
   :maxdepth: 1
   :caption: Introduction
   :hidden:

   intro
   rate_types
   install
   contributing

.. toctree::
   :maxdepth: 1
   :caption: Nuclear Properties
   :hidden:

   sources
   nucleus

.. toctree::
   :maxdepth: 1
   :caption: Exploring Networks in Python
   :hidden:

   pynucastro-examples.ipynb
   plot-types.ipynb
   pynucastro-integration.ipynb


.. toctree::
   :maxdepth: 1
   :caption: Working with Libraries
   :hidden:

   library-examples.ipynb
   tabulated-weak-rates.ipynb

.. toctree::
   :maxdepth: 1
   :caption: Building Networks
   :hidden:

   duplicate-rates.ipynb
   validate-example.ipynb
   electron-capture-example.ipynb
   electron-captures.ipynb
   alternate-rates.ipynb

.. toctree::
   :maxdepth: 1
   :caption: Advanced Rate Operations
   :hidden:

   screening-examples
   modify-example.ipynb
   custom-rates.ipynb
   partition-function

.. toctree::
   :maxdepth: 1
   :caption: Advanced Network Operations
   :hidden:

   approximate-rates
   nse-protons.ipynb
   unimportant-rates.ipynb
   network-cycles.ipynb

.. toctree::
   :maxdepth: 1
   :caption: Nuclear Statistical Equilibrium
   :hidden:

   NSE-example
   nse_table

.. toctree::
   :maxdepth: 1
   :caption: Plasma Neutrinos
   :hidden:

   neutrino-cooling.ipynb

.. toctree::
   :maxdepth: 1
   :caption: Exporting Networks
   :hidden:

   networks
   python-network
   simple-cxx-network
   fortran-network
   amrex-astro-cxx-networks
   julia

.. toctree::
   :maxdepth: 1
   :caption: Some Useful Networks
   :hidden:

   he-burning-example

.. toctree::
   :maxdepth: 1
   :caption: Examples in Nuclear Astrophysics
   :hidden:

   examples/binding-energy.ipynb
   examples/pp-cno.ipynb
   examples/hot-CNO-breakout-example.ipynb
   examples/triple_alpha_eval.ipynb
   examples/he-burning.ipynb
   examples/supernova-lightcurve.ipynb

.. toctree::
   :maxdepth: 1
   :caption: Reference
   :hidden:

   citing
   changes
   API <modules>
   zreferences

.. toctree::
   :maxdepth: 1
   :caption: Index
   :hidden:

   genindex
   modindex
