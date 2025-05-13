Approximate and Modified Rates
==============================

For multidimensional simulation codes, it is common to approximate
some rate sequences to reduce the number of nuclei while keeping the
effective flows.  For example, a common approximation is to combine
$A(\alpha,\gamma)B$ and $A(\alpha,p)X(p,\gamma)$ into an effective
rate and not carry nuclei $X$ in the network.

pynucastro can do several rate approximations.

.. toctree::
   :maxdepth: 1

   alpha-gamma-approx.ipynb
   double-n-capture.ipynb
   modified-rates.ipynb
