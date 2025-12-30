NSE Intro
=========

Nuclear statistical equilibrium (NSE) describes a state in which all strong
nuclear reactions are in chemical equilibrium, with each forward reaction
balanced by its reverse. NSE is achieved when nuclear reaction timescales
are much shorter than hydrodynamic timescales, typically at temperatures
:math:`T \gtrsim 6 \times 10^9\,\mathrm{K}`.

Mathematical Description of NSE
-------------------------------

Imagine that every isotope is connected to the free nucleons,
neutrons and protons, through a series of strong nuclear reactions.
Now since in NSE, every strong reaction is in chemical equilibrium, this means
that the chemical potential of any isotope can be expressed as the sum of the
chemical potentials of the its constituent neutrons and protons.
This can be viewed as a net reaction of the form:

.. math:: (A, Z) \rightleftharpoons Z \  \mathrm{p} + N \ \mathrm{n}

And chemical equilibrium of such reaction is:

.. math:: \mu_i = Z_i \mu_p + N_i \mu_n

Now recall that for gas that follows Maxwell-Boltzmann statistics,
as derived in :ref:`gas_model`, the chemical potential of any isotope is:

.. math:: \mu_i = \mu_i^{\mathrm{kin}} + m_i c^2 + \mu_i^c

Combining with the NSE chemical equilibrium condition, we obtain:

.. math::
   \mu_i^{\mathrm{kin}} = Z_i (\mu_p^{\mathrm{kin}} + \mu_p^c) + N_i \mu_n^{\mathrm{kin}} + \underbrace{(Z_i m_p + N_i m_n - m_i) c^2}_{B_i} - \mu^c_i

Where :math:`B_i` is the binding energy of the isotope. The above expression
can then be used in the mass fraction equation derived from Maxwell-Boltzmann statistics
shown in :ref:`gas_model`, giving the mass abundance equation for system in NSE.

.. math::

   X_{i, \mathrm{NSE}} = \dfrac{(m_i)^{5/2}}{\rho} g_i \left( \dfrac{k_B T}{2\pi \hbar^2} \right)^{3/2} \exp{\left( \dfrac{Z_i (\mu_p^{\mathrm{kin}} + \mu_p^c) + N_i \mu_n^{\mathrm{kin}} + B_i - \mu^c_i}{k_B T} \right)}

NSE Constraint Equations
------------------------

Given the input, :math:`(\rho, T, Y_e)`, the NSE abundance equation has two unknowns,
:math:`\mu_p^{\mathrm{kin}}` and :math:`\mu_n^{\mathrm{kin}}`, which can be solved
subject two constraint equations:

1. Conservation of mass, or a constraint on the input density, :math:`\rho`.

.. math:: \sum_i X_i^{\mathrm{NSE}}(\mu_p^{\mathrm{kin}}, \mu_n^{\mathrm{kin}}) - 1 = 0

2. Conservation of charge for strong reactions, or a constraint on the
   input electron fraction, :math:`Y_e`

.. math:: \sum_i \frac{Z_i}{A_i} X_i^{\mathrm{NSE}} (\mu_p^{\mathrm{kin}}, \mu_n^{\mathrm{kin}}) - Y_e = 0

The corresponding constraint Jacobian is then:

.. math::

   \begin{bmatrix}
   \frac{\partial}{\partial \mu_p^{\mathrm{kin}}} \left(\sum_i X_i^{\mathrm{NSE}} - 1\right) = \sum_i \frac{Z_i}{k_B T} X_i^{\mathrm{NSE}} & \frac{\partial}{\partial \mu_n^{\mathrm{kin}}} \left(\sum_i X_i^{\mathrm{NSE}} - 1\right) = \sum_i \frac{N_i}{k_B T} X_i^{\mathrm{NSE}} \\
   \frac{\partial}{\partial \mu_p^{\mathrm{kin}}} \left(\sum_i \frac{Z_i}{A_i} X_i^{\mathrm{NSE}} - Y_e\right) = \sum_i \frac{Z_i^2}{A_i k_B T} X_i^{\mathrm{NSE}} & \frac{\partial}{\partial \mu_n^{\mathrm{kin}}} \left(\sum_i \frac{Z_i}{A_i} X_i^{\mathrm{NSE}} - Y_e\right) = \sum_i \frac{Z_i N_i}{A_i k_B T} X_i^{\mathrm{NSE}}
   \end{bmatrix}

The NSE abundance system defined by the constraint equations and Jacobian can
be solved using :class:`NSENetwork <pynucastro.networks.nse_network.NSENetwork>`.
The primary interface for this calculation is
:func:`get_comp_nse <pynucastro.networks.nse_network.NSENetwork.get_comp_nse>`.
