Reverse Rates
=============

In this section we will implement the partition functions discussed in
:cite:t:`rauscher:2000` and :cite:t:`rauscher:2003`
with the purpose to compute inverse rates at extreme temperatures and heavy nuclide conditions.

Semi-relativistic ideal gas model
---------------------------------

In order to illustrate our implementation, we will start first by computing the number density of reaction specie in terms of its chemical potential. Let us start our discussion by introducing the 1-particle *canonical partition function* by

.. math:: Z_{1-p} = \dfrac{g}{(2\pi \hbar)^3} \int \, d^3q \,d^3p \, \exp(-E/(kT))

Our purpose is to compute the number density in terms of the chemical potential by assuming a first order velocity approximation of the energy
for each particle as:

.. math:: E = \dfrac{p^2}{2m} + mc^2

where:

.. math:: 
   :nowrap:
   
   \begin{eqnarray}
   Z_{1-p} &=& e^{-mc^2/(kT)}\times\dfrac{g}{(2\pi \hbar)^3} \int \, d^3q \,d^3p \, \exp(-p^2/(2mkT)) \\
           &=& e^{-mc^2/(kT)}\times gV\left(\dfrac{kT}{2\pi\hbar^2} \right)^{3/2}
   \end{eqnarray}

From this point, we may compute the N-particle *canonical partition function* by

.. math:: Z_N = \dfrac{1}{N!}Z_{1-p}^N

and by using the fugacity parameter, defined by :math:`z=e^{\mu/(kT)}`, we also may compute the *grand-canonical partition function*:

.. math:: \mathcal{Z} = \sum_{N=0} z^N Z_N = \sum_{N=0} \dfrac{(zZ_{1-p})^N}{N!} = \exp\left[e^{-mc^2/(kT)}\times zgV\left(\dfrac{kT}{2\pi \hbar^2} \right)^{3/2} \right]

with the associated grand-canonical potential:

.. math:: \Omega = -kT\log \mathcal{Z} = -kT gV\left(\dfrac{kT}{2\pi \hbar^2} \right)^{3/2} e^{\mu/(kT)}e^{-mc^2/(kT)}

and the number of particles may be computed as the first derivative of the grand-canonical potential with respect to the chemical potential :math:`\mu`:

.. math:: N = - \left( \dfrac{\partial \Omega}{\partial \mu} \right ) = gVm^{3/2}\left( \dfrac{kT}{2\pi \hbar^2} \right)^{3/2} e^{\mu/(kT)}e^{-mc^2/(kT)}

From the number density definition, :math:`n=N/V` we may compute finally, the relationship between the number density and the chemical potential for a
semi-relativistic ideal gas by:

.. math:: n = gm^{3/2} \left( \dfrac{kT}{2\pi \hbar^2} \right)^{3/2} e^{\mu/(kT)}e^{-mc^2/(kT)}

or equivalently,

.. math:: \mu = kT \log \left[\dfrac{n}{g}\left( \dfrac{2\pi \hbar^2}{mkT} \right)^{3/2} \right] + mc^2

as summarized in :cite:`rauscher:2000`.

Reverse rate definition
-----------------------

The forward and reverse reactions are entangled due to thermodynamical and nuclear equilibrium aspects. In this section we will see how this connection is introduced
and explore its consequences in the nuclear rates computation, as discussed in :cite:`rauscher:2000`. Starting with the reaction:

.. math:: A + B \rightarrow C

characterized by a positive Q-capture energy, we may define the reverse reaction as:

.. math:: C \rightarrow A + B

clearly characterized by a negative Q-capture energy. (It releases energy after the decay of C). Therefore, the rate of A and B particles in the forward reaction implies:

.. math:: \dfrac{dn_{AB}}{dt} = -n_An_B\langle \sigma v \rangle_{AB}

on the reverse reaction, we have:

.. math:: \dfrac{dn_C}{dt} = -\lambda_C n_C

Now, taking in consideration that in equilibrium, two equations are satisfied. First, the conservation of nucleons and second, the chemical potential balance:

.. math::

   \dfrac{dn_{AB}}{dt} - \dfrac{dn_c}{dt} = 0



   \mu_A + \mu_B = \mu_C

From, the first equation we get the following relationship between :math:`\lambda_C` and :math:`\langle \sigma v \rangle_{AB}` :

.. math:: \dfrac{n_A n_B}{n_C} = \dfrac{\lambda_C}{\langle \sigma v \rangle_{AB}}

Now, the left term of the previous equation may be computed by using the definition computed in the previous section, and the second equation:

.. math::
   :nowrap:

   \begin{eqnarray}
   \dfrac{n_An_B}{n_C}&=&\left(\dfrac{m_Am_B}{m_C} \right)^{3/2}\dfrac{g_Ag_B}{g_C} \left(\dfrac{kT}{2\pi \hbar^2} \right)^{3/2} \times e^{(\mu_A+\mu_B-\mu_C)/(kT)} \times e^{-(m_A+m_B-m_C)c^2/(kT)}\\
                      &=&\left(\dfrac{m_Am_B}{m_C} \right)^{3/2}\dfrac{g_Ag_B}{g_C} \left(\dfrac{kT}{2\pi \hbar^2} \right)^{3/2} \times e^{-Q/(kT)}\\
                      &=&\left(\dfrac{A_AA_B}{A_C} \right)^{3/2}\dfrac{g_Ag_B}{g_C} \left(\dfrac{m_ukT}{2\pi \hbar^2}  \right)^{3/2} \times e^{-Q/(kT)}
   \end{eqnarray}

Therefore, we can compute the photodisintegration rates using the capture rates by:

.. math:: 
   :nowrap:

   \begin{eqnarray}
   \dfrac{\lambda_C}{N_a\langle \sigma v \rangle_{AB}} &= \left(\dfrac{A_AA_B}{A_C} \right)^{3/2}\dfrac{g_Ag_B}{g_C} \left(\dfrac{m_ukT}{2\pi \hbar^2}  \right)^{3/2}\dfrac{1}{N_a} \times e^{-Q/(kT)}\\
                                                       &= \left(\dfrac{A_AA_B}{A_C} \right)^{3/2}\dfrac{g_Ag_B}{g_C} T^{3/2}F \times e^{-Q/(kT)}
   \end{eqnarray}

where the numerical factor :math:`F` is defined by:

.. math:: T^{3/2}F = \left(\dfrac{m_ukT}{2\pi \hbar^2}  \right)^{3/2}\dfrac{1}{N_a}

Similarly, for a forward reaction:

.. math:: A + B \rightarrow C + D

we can compute the reverse rate :math:`\langle \sigma v \rangle_{CD}` in terms of the forward rate :math:`\langle \sigma v \rangle_{AB}` as

.. math:: \dfrac{N_a\langle \sigma v \rangle_{CD}}{N_a\langle \sigma v \rangle_{AB}} = \left(\dfrac{A_AA_B}{A_CA_D} \right)^{3/2}\dfrac{g_Ag_B}{g_Cg_D}  \times e^{-Q/(kT)}


Partition functions
-------------------

Until now, we have not discussed the role of :math:`g`, which encompasses the number of spin states that a particle may adopt. In reactions :math:`i(j,o)m`  the target :math:`i` and residual nucleus :math:`m` contributes significantly more in the calculation of the rates than the incident and outgoing particles due to their complexity and the number of quantum levels they may assume. Therefore we may consider :math:`\langle \sigma v \rangle_{ij} \rightarrow \langle \sigma v \rangle_i`  and  :math:`\langle \sigma v \rangle_{om} \rightarrow \langle \sigma v \rangle_m`, which symbolize the nucleus :math:`i`, or :math:`m`, and all the remaining particles in their channels. 

Using this notation, capture reactions rates of the type :math:`(p,\gamma)`, :math:`(n,\gamma)`, :math:`(\alpha,\gamma)`, and their reverse photodisintegration decay rates, are then related by:

.. math:: \dfrac{\lambda_{\gamma}}{N_a\langle \sigma v \rangle_i} = \left(\dfrac{A_iA_j}{A_m} \right)^{3/2}\dfrac{g_ig_j}{g_m} T^{3/2}F \times e^{-Q/(kT)}

otherwise, the forward and reaction rates are related by:

.. math:: \dfrac{N_a\langle \sigma v \rangle_m}{N_a\langle \sigma v \rangle_i} = \left(\dfrac{A_iA_j}{A_oA_m} \right)^{3/2}\dfrac{g_ig_j}{g_og_m}  \times e^{-Q/(kT)}

Now, the ingoing and outgoing are far apart in comparison with the nucleus size before the collision takes place; this allow us to assume that particles like :math:`j` and :math:`o` are in their ground state as they approach. Hence,
:math:`g_j = (2J_j+1)` and :math:`g_o = (2J_o + 1)` where :math:`J_j` and :math:`J_o` are the :math:`j` and :math:`o` particle spin in their respective ground state. However, the target and residual nucleus :math:`i` and :math:`m` particle spin, may be in a superposition of states due to excitation of upper levels caused by an steady increase of the temperature. As pointed out in :cite:`rauscher:2000`, we have to replace :math:`g_i\rightarrow (2J_i+1)G_i` and :math:`g_m\rightarrow (2J_m+1)G_m`, where

.. math::
   :nowrap:

   \begin{eqnarray}
   (2J_i+1)G_i &=& \sum_{\mu} (2J^{\mu}_i+1) \exp\left(-E_i^{\mu}/(kT)\right)\\
   (2J_m+1)G_m &=& \sum_{\nu} (2J^{\nu}_m+1) \exp\left(-E_m^{\nu}/(kT)\right)
   \end{eqnarray}

The quantities :math:`G_i` and :math:`G_j` are the target and residual partition functions that are normalized with respect to their ground state particle spin :math:`J_i` and :math:`J_m` respectively. Now, we are in position to write the following relationships between the forward and reverse rates:

.. math::
   :nowrap:

   \begin{eqnarray}
   \lambda_{\gamma}{N_A \langle \sigma v \rangle_i} &=& \left(\dfrac{A_iA_j}{A_m} \right)^{3/2}\dfrac{(2J_i+1)(2J_j+1)}{(2J_m+1)} \dfrac{G_i}{G_m} T^{3/2}F \times e^{-Q/(kT)}\\
   \dfrac{N_a\langle \sigma v \rangle_m}{N_a\langle \sigma v \rangle_i} &=& \left(\dfrac{A_iA_j}{A_oA_m} \right)^{3/2}\dfrac{(2J_i+1)(2J_j+1)}{(2J_o+1)(2J_m+1)} \dfrac{G_i}{G_m}  \times e^{-Q/(kT)}
   \end{eqnarray}
   
or equivalently, after absorbing all the quantities in the forward rate with the exception of :math:`G_i`:

.. math::
   :nowrap:
   
   \begin{eqnarray}
   \lambda_{\gamma} &=& \lambda_{\gamma}'\dfrac{G_i}{G_m}\\
   N_a\langle \sigma v \rangle_m &=& N_a\langle \sigma v \rangle_m' \dfrac{G_i}{G_m}
   \end{eqnarray}
   
where the quantities :math:`\lambda_{\gamma}'` and :math:`N_a\langle \sigma v \rangle_m'` are provided by REACLIB, under the ``-v`` flag.
   
Implementing partition functions
--------------------------------

The partition function information is contained in three main classes:

* :class:`PartitionFunction <pynucastro.nucdata.partition_function.PartitionFunction>` materialize the temperature and the partition function values into an object, interpolating across all the defined points using a cubic spline interpolation. If a temperature value is outside the temperature range, we keep its value constant to the nearest boundary value. 

* :class:`PartitionFunctionTable <pynucastro.nucdata.partition_function.PartitionFunctionTable>`  reads a table and construct a dictionary between each nucleus and their partition function class object.

* :class:`PartitionFunctionCollection <pynucastro.nucdata.partition_function.PartitionFunctionCollection>` collects all the formatted table information, inside ``/nucdata/PartitionFunction/``. It allow us to include the high temperature tables in :cite:`rauscher:2003` and to select the model used to compute the partition functions, respectively. By default, we include high temperatures, and our partition function model to be the *finite range droplet model (FRDM)*. If a nucleus is not in the collection, we set the partition function values to 1.0 by default.

Inside the :class:`Nucleus <pynucastro.nucdata.nucleus.Nucleus>` class, we have defined ``set_partition_function()`` which setup our partition function collection, our high temperatures consideration, and the model used to compute the partition function data. On the other hand, ``get_partition_function()`` assigns a partition function class object to the defined nucleus. Let us illustrate how it work:

.. code-block:: python
   
   import pynucastro
   
   co46 = pynucastro.rates.Nucleus('co46')
   pCollection = pynucastro.nucdata.PartitionFunctionCollection()

   co46.set_partition_function(pCollection=pCollection, set_data='etfsiq', use_high_temperatures=True)
   pf_co46 = co46.get_partition_function()

Now, from this point we define a method inside a :class:`Rate <pynucastro.rates.rate.Rate>` named ``set_partition_function()`` which reads a partition function collection
and setup all the nucleus inside the reaction rate. Let us illustrate now, how it works:

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
 



   

   

 
