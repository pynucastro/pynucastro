pynucastro.networks.base\_cxx\_network module
=============================================

The generated C++ Jacobian includes explicit composition derivatives of
rate coefficients for species listed in ``Rate.rate_comp_dependence``.
For a reaction contribution :math:`f_j = C_j(Y, \rho)\lambda(Y)`, it uses
the product rule,
:math:`\partial f_j/\partial Y_i = \lambda\,\partial C_j/\partial Y_i
+ C_j\,\partial\lambda/\partial Y_i`.
The coefficient derivatives are read from ``rate_derivs_t``; for example,
the effective Fe52(nn, gamma)Fe54 rate contributes a term using
``drate_Fe52_n_n_to_Fe54_approx_dYN``.

.. automodule:: pynucastro.networks.base_cxx_network
   :members:
   :undoc-members:
   :show-inheritance:
