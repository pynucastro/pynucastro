Derived Rates
=============

In this section we discuss the implementation of
:class:`DerivedRate <pynucastro.rates.derived_rate.DerivedRate>`
via detailed balance.
We will also discuss the partition functions mentioned in
:cite:t:`rauscher:2000` and :cite:t:`rauscher:2003`
with the purpose to compute inverse rates at extreme temperatures and
heavy nuclide conditions.

.. _gas_model:

Semi-relativistic ideal gas model
---------------------------------

In order to illustrate our implementation, we will start first by computing
the number density of reaction specie in terms of its chemical potential.
Let us start our discussion by introducing the 1-particle *canonical partition function* by

.. math:: Z_{1-p} = \sum_l (2J_l + 1) \dfrac{1}{(2\pi \hbar)^3} \int \, d^3q \,d^3p \, \exp{\left(-\dfrac{\epsilon(l, \mathbf{p}, \mathbf{q})}{k_B T}\right)}

Our purpose is to compute the number density in terms of the chemical potential
by assuming a first order velocity approximation of the energy, internal excitation energy,
rest energy and coulomb screening potential for each particle as:

.. math:: \epsilon(l, \mathbf{p}, \mathbf{q}) = \underbrace{\dfrac{p^2}{2m}}_{\mathrm{Kinetic}} + \underbrace{\Delta_l}_{\mathrm{Excitation}} + \underbrace{mc^2}_{\mathrm{Rest}} + \underbrace{\mu^c}_{\mathrm{Coulomb}}

.. math::
   :nowrap:

   \begin{eqnarray}
   Z_{1-p} &=& \underbrace{\sum_l(2J_l + 1) e^{-\Delta_l/(k_B T)}}_{g} \dfrac{1}{(2\pi \hbar)^3} \int \, d^3q \,d^3p \, \exp{\left(-\dfrac{\mathbf{p}^2/2m + mc^2 + \mu^c}{k_B T} \right)} \\
           &=& e^{-(mc^2 + \mu^c)/(k_B T)} \times gV\left(\dfrac{m k_B T}{2\pi\hbar^2} \right)^{3/2}
   \end{eqnarray}

From this point, we may compute the N-particle *canonical partition function* by

.. math:: Z_N = \dfrac{1}{N!}Z_{1-p}^N

and by using the fugacity parameter, defined by :math:`z=e^{\mu/(kT)}`,
we also may compute the *grand-canonical partition function*:

.. math:: \mathcal{Z} = \sum_{N=0}^{\infty} z^N Z_N = \sum_{N=0}^{\infty} \dfrac{(zZ_{1-p})^N}{N!} = \exp\left[e^{-(mc^2+\mu^c)/(kT)}\times zgV\left(\dfrac{k_B T}{2\pi \hbar^2} \right)^{3/2} \right]

with the associated grand-canonical potential:

.. math:: \Omega = -k_B T\log \mathcal{Z} = -k_B T gV\left(\dfrac{k_B T}{2\pi \hbar^2} \right)^{3/2} e^{(\mu - mc^2 - \mu^c)/(k_B T)}

and the number of particles may be computed as the first derivative of the
grand-canonical potential with respect to the chemical potential :math:`\mu`:

.. math:: N = - \left( \dfrac{\partial \Omega}{\partial \mu} \right ) = gV\left( \dfrac{m k_B T}{2\pi \hbar^2} \right)^{3/2} e^{(\mu - mc^2 - \mu^c)/(k_B T)}

From the number density definition, :math:`n=N/V` we may compute finally,
the relationship between the number density and the chemical potential for a
semi-relativistic ideal gas by:

.. math:: n = g \left( \dfrac{m k_B T}{2\pi \hbar^2} \right)^{3/2} e^{(\mu - mc^2 - \mu^c)/(k_B T)}

or equivalently,

.. math:: \mu = \underbrace{k_B T \log \left[\dfrac{n}{g}\left( \dfrac{2\pi \hbar^2}{mkT} \right)^{3/2} \right]}_{\mu^{kin}} + mc^2 + \mu^c

as summarized in :cite:`rauscher:2000`.

Often times we want to work in massfractions, :math:`X_i = m_i n_i / \rho`,
or in molar fractions, :math:`Y_i = X_i / A_i`. Hence:

.. math:: Y_i = \dfrac{1}{\rho} \dfrac{(m_i)^{5/2}}{A_i} g_i \left( \dfrac{k_B T}{2\pi \hbar^2} \right)^{3/2} e^{(\mu_i - m_i c^2 - \mu_i^c)/(k_B T)}

Partition functions
-------------------

Until now, we have not discussed the role of :math:`g`, which encompasses the
number of spin states that a particle may adopt. Now, since the ingoing and
outgoing particles are far apart in comparison with the nucleus size before
the collision takes place; this allow us to assume that particles are likely
in their ground state as they approach. Hence to account for state degeneracy,
:math:`g = 2J_0 + 1`, where :math:`J_0` is the particle spin in the ground state.
However, particles may be in a superposition of states due to excitation
of upper levels, with excitation energies :math:`\Delta_l`,
caused by an steady increase of the temperature.
As pointed out in :cite:`rauscher:2000`, we have to replace
:math:`g \rightarrow (2J_0 + 1) G`, where :math:`G` is the internal partition
function normalized by the ground state spin.

.. math::
   g \equiv (2 J_0 + 1) \times \underbrace{\frac{\sum_l(2J_l + 1) e^{-\Delta_l/(k_B T)}}{2 J_0 + 1}}_{G}

.. note::

   ReacLib rates with ``-v`` flag are derived similarly, but setting :math:`G = 1`.

Implementing partition functions
--------------------------------

The partition function information is contained in three main classes:

* :class:`PartitionFunction <pynucastro.nucdata.partition_function.PartitionFunction>`
  materialize the temperature and the partition function values into an object,
  interpolating across all the defined points using a cubic spline interpolation.
  If a temperature value is outside the temperature range, we keep its value
  constant to the nearest boundary value.

* :class:`PartitionFunctionTable <pynucastro.nucdata.partition_function.PartitionFunctionTable>`
  reads a table and construct a dictionary between each nucleus and their partition function class object.

* :class:`PartitionFunctionCollection <pynucastro.nucdata.partition_function.PartitionFunctionCollection>`
  collects all the formatted table information, inside ``/nucdata/PartitionFunction/``.
  It allow us to include the high temperature tables in :cite:`rauscher:2003` and to
  select the model used to compute the partition functions, respectively. By default,
  we include high temperatures, and our partition function model to be the
  *finite range droplet model (FRDM)*. If a nucleus is not in the collection,
  we set the partition function values to 1.0 by default.

Inside the :class:`Nucleus <pynucastro.nucdata.nucleus.Nucleus>` class,
we have defined ``set_partition_function()`` which setup our partition function
collection, our high temperatures consideration, and the model used to compute
the partition function data. On the other hand, ``get_partition_function()``
assigns a partition function class object to the defined nucleus.
Let us illustrate how it work:

.. code-block:: python

   import pynucastro

   co46 = pynucastro.rates.Nucleus('co46')
   pCollection = pynucastro.nucdata.PartitionFunctionCollection()

   co46.set_partition_function(pCollection=pCollection, set_data='etfsiq', use_high_temperatures=True)
   pf_co46 = co46.get_partition_function()

Now, from this point we define a method inside a :class:`Rate <pynucastro.rates.rate.Rate>`
named ``set_partition_function()`` which reads a partition function collection
and setup all the nucleus inside the reaction rate. Let us illustrate now,
how it works:

.. code-block:: python

   import pynucastro

   p_collection = pynucastro.nucdata.PartitionFunctionCollection()

   o18_pg_f19 = pynucastro.rates.Rate('../library/o18-pg-f19-il10')
   o18_pg_f19.set_partition_function(p_collection=p_collection, set_data='frdm', use_high_temperatures=True)

   p = o18_pg_f19.reactants[0]
   o18 = o18_pg_f19.reactants[1]
   f19 = o18_pg_f19.products[0]

   pf_p = p.get_partition_function()
   pf_o18 = o18.get_partition_function()
   pf_f19 = f19.get_partition_function()


Derived rate definition
-----------------------

The forward and reverse reactions are entangled due to thermodynamical and
nuclear equilibrium aspects. In this section we will see how this connection is introduced
and explore its consequences in the nuclear rates computation,
as discussed in :cite:`rauscher:2000`. The general case is worked out in the appendix of
:cite:`smith_pynucastro_2023`. Here we demonstrate a less general reaction channel,
:math:`A + B \rightarrow C + D`, discussed in :cite:`ReacLib`, for simplicity.

Starting with a pair of reaction channels, where $A$, $B$, $C$, $D$ are distinct
nuclear species:

.. math:: N_A \  A + N_B \ B \rightleftharpoons N_C \ C + N_D \ D

The coupled ordinary differential equations characterizing the evolution
of the species can be written as the following:

.. math::

    \begin{eqnarray}
      \frac{1}{N_A} \frac{d}{d t} Y_A &= \frac{1}{N_B} \frac{d}{d t} Y_B &= - \frac{\rho^{\nu_r} Y_A^{N_A} Y_B^{N_B}}{N_A! N_B!} \underbrace{N_a^{\nu_r} \langle \sigma v \rangle}_{\lambda_\mathrm{for}} + \frac{\rho^{\nu_p} Y_C^{N_C} Y_D^{N_D}}{N_C! N_D!} \underbrace{N_a^{\nu_p} \langle \sigma v \rangle}_{\lambda_\mathrm{rev}} \\
      \frac{1}{N_C} \frac{d}{d t} Y_C &= \frac{1}{N_D} \frac{d}{d t} Y_D &= + \frac{\rho^{\nu_r} Y_A^{N_A} Y_B^{N_B}}{N_A! N_B!} \lambda_\mathrm{for} - \frac{\rho^{\nu_p} Y_C^{N_C} Y_D^{N_D}}{N_C! N_D!} \lambda_\mathrm{rev}
    \end{eqnarray}

where $\nu_r = N_A + N_B - 1$ and $\nu_p = N_C + N_D - 1$
Now, taking in consideration that forward and reverse rates are in equilibrium:

.. math::

   \frac{d }{d t} Y_A = \frac{d }{d t} Y_B = \frac{d }{d t} Y_C = \frac{d }{d t} Y_D = 0

Hence the ratio between the screend forward and the screened reverse rate can be written as:

.. math::

   \frac{\lambda^\mathrm{scn}_\mathrm{rev}}{\lambda^\mathrm{scn}_\mathrm{for}} = \rho^{\nu_r - \nu_p} \ \frac{Y^{N_A}_{\mathrm{eq},A} Y^{N_B}_{\mathrm{eq},B}}{Y^{N_C}_{\mathrm{eq},C} Y^{N_D}_{\mathrm{eq},D}} \ \frac{N_C! N_D!}{N_A! N_B!}

To continue further, the molar fraction at equilibrium state is characterized by
the form of molar fraction that we obtained in :ref:`gas_model` with
the chemical potential representing the equilibrium state,
i.e. :math:`\mu = \mu^{\mathrm{eq}}`. With the definition of the Q-value of the
reverse rate:

.. math::

   Q_\mathrm{rev} \equiv (N_C \ m_C + N_D \ m_D) c^2 - (N_A \ m_A + N_B \ m_B) c^2

And the detailed balance condition, i.e. chemical equilibrium:

.. math::

     N_A \ \mu^\mathrm{eq}_A + N_B \ \mu^\mathrm{eq}_B = N_C \ \mu^\mathrm{eq}_C + N_D \ \mu^\mathrm{eq}_D

We finally obtain:

.. math::

   \frac{\lambda^\mathrm{scn}_\mathrm{rev}}{\lambda^\mathrm{scn}_\mathrm{for}} = \underbrace{\frac{K_A^{N_A} K_B^{N_B}}{K_C^{N_C} K_D^{N_D}} \frac{N_C! N_D!}{N_A! N_B!} \exp{\left(\frac{Q_{\mathrm{rev}}}{k_B T}\right)}}_{R} \exp{ \left(-\frac{N_A \ \mu^c_A + N_B \ \mu^c_B - N_C \ \mu^c_C - N_D \ \mu^c_D}{k_B T} \right)}

where :math:`R` is the equilibrium ratio without the screening term
and the prefactor :math:`K` is defined as:

.. math::

   K_i \equiv \frac{m_i^{5/2}}{A_i} g_i \left(\frac{k_B T}{2\pi\hbar^2}\right)^{3/2}

Now assuming that the electron screening routine follows linear mixing rule
:cite:`dewitt_screening_1973`, then the screening enhancement factors for
the forward and reverse reaction channels that we're considering where
:math:`\lambda^\mathrm{scn}_\mathrm{rev} = \lambda^\mathrm{uns}_\mathrm{rev} \ F^\mathrm{scn}_\mathrm{rev}`
and :math:`\lambda^\mathrm{scn}_\mathrm{for} = \lambda^\mathrm{uns}_\mathrm{for} \ F^\mathrm{scn}_\mathrm{for}`
will be:

.. note::

   The only screening routine that follows the linear mixing rule is
   :func:`potekhin_1998 <pynucastro.screening.screen.potekhin_1998>`.
   This is important if you want to be self-consistent with NSE.


.. math::

   \begin{eqnarray}
   F^\mathrm{scn}_\mathrm{for} &= \exp{\left(\frac{N_A \mu^c(Z_A) + N_B \mu^c(Z_B) - \mu^c(N_A Z_A + N_B Z_B)}{k_B T}\right)} \\
   F^\mathrm{scn}_\mathrm{rev} &= \exp{\left(\frac{N_C \mu^c(Z_C) + N_D \mu^c(Z_D) - \mu^c(N_C Z_C + N_D Z_D)}{k_B T}\right)}
   \end{eqnarray}

where :math:`Z` is the atomic number of the reactants involved in the reaction.
With the above rule, the ratio between the screened reverse rate and the
unscreened forward rate is:

.. math::

   \frac{\lambda^\mathrm{scn}_\mathrm{rev}}{\lambda^\mathrm{uns}_\mathrm{for}} = R \times \exp{ \left(\frac{N_C \ \mu^c_C + N_D \ \mu^c_D - \mu^c(N_A Z_A + N_B Z_B) }{k_B T} \right)}

Now due to conservation of charge for strong reactions:

.. math::

   N_A Z_A + N_B Z_B = N_C Z_C + N_D Z_D

We finally obtain:

.. math::

   \frac{\lambda^\mathrm{scn}_\mathrm{rev}}{\lambda^\mathrm{uns}_\mathrm{for}} = R \ F^\mathrm{scn}_\mathrm{rev}

We note the above expression is self-consistent with the definition
of the electron screening enhancement where:

.. math::

   \lambda^\mathrm{scn}_\mathrm{rev} = \underbrace{\lambda^\mathrm{uns}_\mathrm{for} \ R}_{\lambda^\mathrm{uns}_\mathrm{rev}} \ F^\mathrm{scn}_\mathrm{rev}


Therefore, the unscreened
:class:`DerivedRate <pynucastro.rates.derived_rate.DerivedRate>`
is evaluated by taking the provided source rate and multiplying it by the
temperature-dependent equilibrium ratio, :math:`R`. The appropriate
screening enhancement factor can then be computed and applied in the usual way,
just as for any other rate.

.. note::

   It is completely general to compute the derived *forward* rate
   based on the reverse rate. Examples of this are forward rates,
   i.e. ``Q > 0``, provided by ReacLib but with the ``-v`` flag.
