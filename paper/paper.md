---
title: 'pynucastro: an interface to nuclear reaction rates and code generator for reaction network equations'

tags:
- nuclear
- reactions
- astrophysics
- physics
- visualization
- code generation

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
at a particular temperature, and find the exponent, n, for a simply
$\rm{T^n}$ parameterization.  From a collection of rates, the flow between
the nuclei can be visualized interactively using Jupyter widgets.
These features help both with designing a network for a simulation as
well as for teaching nuclear astrophysics in the classroom.

After selecting a set of rates for a given problem, pynucastro can
construct a reaction network from those rates, using SymPy
[@SymPy.2017] to determine the system of ODEs comprising the
network. From the symbolic expressions for the ODE right hand side,
pynucastro can then generate Python or Fortran code to calculate the
right hand side. Generated Python right hand sides evolve species in
the reaction network and pynucastro includes a Python example
integrating the CNO cycle for hydrogen burning.

For Fortran networks, pynucastro will also generate a routine to
compute the analytic jacobian matrix for implicit integration. These
networks incorporate weak, intermediate, and strong reaction rate
screening for the Reaclib rates [@Screening.Graboske.1973;
@Screening.Alastuey.1978; @Screening.Itoh.1979]. Fortran networks can
also incorporate selected weak reaction rate tabulations
[@Suzuki.2016]. To calculate energy generation in Fortran networks,
pynucastro incorporates nuclear binding energies from the Atomic Mass
Data Center [@AME2016.1; @AME2016.2] and the 2014 CODATA recommended
values for the fundamental physical constants [@CODATA.2014].

pynucastro is capable of generating two kinds of Fortran reaction
networks. The first type is a standalone network with a driver program
to integrate species and energy generation using the variable-order
ODE integration package VODE [@VODE.1989]. This Fortran driver program
is designed to be easy to use and can integrate reaction networks
significantly faster than is possible for the generated Python
networks.

Secondly, pynucastro can generate a Fortran network consisting of
right hand side and jacobian modules that evolve species, temperature,
and energy generation for the StarKiller Microphysics code. Via
StarKiller Microphysics, astrophysical simulation codes such as Castro
and Maestro can directly use pynucastro reaction networks. pynucastro
includes a carbon-burning network with tabulated $\rm{A=23}$ Urca weak
reactions currently used for research with Maestro
[@Zingale.astronum.2017].

# References
