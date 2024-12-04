Installation
============

PyPI
----

The easiest way to install pynucastro is via pip:

.. prompt:: bash

   pip install pynucastro

.. note::

   If you are using python 3.13, then you will need to install
   the release candidate of numba until the final release is out.
   This can be done as:

   .. prompt:: bash

      pip install numba==0.61.0rc1


Conda
-----

If you use conda, you can install pynucastro from `conda-forge
<https://anaconda.org/conda-forge/pynucastro>`_ as:

.. prompt:: bash

   conda install conda-forge::pynucastro


Building from source
--------------------

Alternately, you can clone the source and install it locally via:

.. prompt:: bash

   git clone https://github.com/pynucastro/pynucastro.git
   cd pynucastro
   pip install --user .
