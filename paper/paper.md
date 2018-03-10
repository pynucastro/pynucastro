---
title: 'pynucastro: an interface to nuclear reaction rates and code generator for reaction network equations'

tags:
- nuclear
- reactions
- astrophysics
- physics
- visualization
- code generation
- integration
- differential equations

authors:
- name: Donald E. Willcox
  orcid: 0000-0003-2300-5165
  affiliation: 1
- name: Michael Zingale
  orcid: 0000-0001-8401-030X
  affiliation: 1

affiliations:
- name: Department of Physics and Astronomy, Stony Brook University
  index: 1

date: 27 January 2018

bibliography: paper.bib
---

# Summary

pynucastro addresses two needs in the field of nuclear astrophysics:
visual exploration of nuclear reaction rates or networks and automated
code generation for integrating reaction network ODEs. pynucastro
accomplishes this by interfacing with nuclear reaction rate
parameterizations published by the JINA Reaclib project
[@Reaclib.2010].

Interactive exploration is enabled by a set of classes that provide
methods to visualize the temperature dependency of a rate, evaluate it
at a particular temperature, and find the exponent, n, for a simple
$\rm{T^n}$ parameterization.  From a collection of rates, the flow
between the nuclei can be visualized interactively using Jupyter
widgets.  These features help both with designing a network for a
simulation as well as for teaching nuclear astrophysics in the
classroom.

After selecting a set of rates for a given problem, pynucastro can
construct a reaction network from those rates consisting of Python
code to calculate the ODE right hand side. Generated Python right hand
sides evolve species in the reaction network, and pynucastro includes
a Python example integrating the CNO cycle for hydrogen burning.

pynucastro can also generate Fortran code implementing reaction
networks, using SymPy [@SymPy.2017] to determine the system of ODEs
comprising the network. From the symbolic expressions for the ODE
right hand side, pynucastro also generates a routine to compute the
analytic Jacobian matrix for implicit integration.

Fortran networks incorporate weak, intermediate, and strong reaction
rate screening for the Reaclib rates [@Screening.Graboske.1973;
@Screening.Alastuey.1978; @Screening.Itoh.1979]. These networks can
also include selected weak reaction rate tabulations
[@Suzuki.2016]. To calculate energy generation in Fortran networks,
pynucastro uses nuclear binding energies from the Atomic Mass Data
Center [@AME2016.1; @AME2016.2] and the 2014 CODATA recommended values
for the fundamental physical constants [@CODATA.2014].

pynucastro is capable of generating two kinds of Fortran reaction
networks. The first type is a standalone network with a driver program
to integrate species and energy generation using the variable-order
ODE integration package VODE [@VODE.1989]. This Fortran driver program
is designed to be easy to use and can integrate reaction networks
significantly faster than is possible for the generated Python
networks.

Secondly, pynucastro can generate a Fortran network consisting of
right hand side and Jacobian modules that evolve species, temperature,
and energy generation for the StarKiller Microphysics code. Via
StarKiller Microphysics, astrophysical simulation codes such as Castro [@castro]
and Maestro [@MAESTRO:Multilevel] can directly use pynucastro reaction networks.
pynucastro includes a carbon burning network with tabulated $\rm{A=23}$ Urca weak
reactions currently used for studying white dwarf convection with
Maestro [@Zingale.astronum.2017].

Future work will focus on implementing nuclear partition functions to
compute reverse reaction rates in the Reaclib library [@RaTh2000;
@Rauscher2003]. It is also in some cases necessary to compute reverse
reaction rates using detailed balance with a consistent nuclear mass
model instead of using the parameterized reverse reaction rates in
Reaclib [@SkyNet.2017]. Additionally, work is ongoing to port the
networks generated for StarKiller Microphysics to CUDA Fortran to
support parallel reaction network integration on GPU systems
[@Zingale.astronum.2017]. We intend to implement this port directly
into the pynucastro-generated networks.

We wish to thank Abigail Bishop for discussions about code generation
for the StarKiller Microphysics code as well as for exploratory
calculations. We are grateful to Max P. Katz for numerous discussions
that enabled the ongoing port of pynucastro-generated networks to CUDA
Fortran. We also wish to thank Christopher Malone for discussions
about various implementation details in pynucastro as well as sample
code to improve element identification. We especially thank Josiah
Schwab for helpful discussions about nuclear partition functions and
reverse rates as well as for testing pynucastro and pointing out
issues in visualization and documentation. This work was supported by
DOE/Office of Nuclear Physics grant DE-FG02-87ER40317.

# References
