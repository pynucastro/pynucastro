"""
Fundamental constants used throughout pynucastro.  Most of these
come from ``scipy.constants``.

"""

import math

import scipy.constants as physical_constants

#: atomic mass unit in g
m_u = physical_constants.value("unified atomic mass unit") * physical_constants.kilo

#: atomic mass unit in MeV
m_u_MeV = physical_constants.value("atomic mass constant energy equivalent in MeV")

#: proton mass in g
m_p = physical_constants.value("proton mass") * physical_constants.kilo

#: proton mass in MeV
m_p_MeV = physical_constants.value("proton mass energy equivalent in MeV")

#: neutron mass in g
m_n = physical_constants.value("neutron mass") * physical_constants.kilo

#: neutron mass in MeV
m_n_MeV = physical_constants.value("neutron mass energy equivalent in MeV")

#: electron mass in g
m_e = physical_constants.value("electron mass") * physical_constants.kilo

#: electron mass in MeV
m_e_MeV = physical_constants.value("electron mass energy equivalent in MeV")

#: Avogadro's number
N_A = physical_constants.value("Avogadro constant")

#: Boltzmann constant in erg/K
k = physical_constants.value("Boltzmann constant") / physical_constants.erg

#: speed of light in cm / s
c_light = physical_constants.value("speed of light in vacuum") / physical_constants.centi

#: Stefan-Boltzmann constant in erg/cm^2/K^4/s
sigma = physical_constants.value("Stefan-Boltzmann constant") / physical_constants.erg * physical_constants.centi**2

#: radiation constant in erg/cm^3/K^4
a = 4.0 * sigma / c_light

# SciPy 1.15.0 fixed how exact values are calculated, which makes them slightly
# more accurate compared to 1.14.1 and earlier. Calculate this one manually for
# consistency between versions (should match 1.15.0).
# This calculation comes directly from scipy.constants._codata, and works
# because 1 eV = 1 V * q_e = 1 J/C * 1.602e-19 C = 1.602e-19 J.

#: Boltzmann constant in MeV/K
k_MeV = physical_constants.value("Boltzmann constant") / physical_constants.e / physical_constants.mega

#: Planck constant in erg s
h = physical_constants.value("Planck constant") / physical_constants.erg

# another exact constant

#: Planck constant / 2 pi in erg s
hbar = h / (2*math.pi)

#: elementary charge in stat-Coulomb
q_e = physical_constants.value("elementary charge") * (physical_constants.c * 100) / 10  # C to statC (esu)

# conversions

#: conversion factor from erg to MeV
erg2MeV = physical_constants.erg / (physical_constants.eV * physical_constants.mega)

#: conversion factor from MeV to erg
MeV2erg = (physical_constants.eV * physical_constants.mega) / physical_constants.erg

# CODATA 2018 values of masses needed for Nubase 2020 consistency
# these should no longer be needed when the next Nubase is released and we drop
# compatibility for SciPy <1.15

#: atomic mass unit in MeV from CODATA 2018
m_u_MeV_C18 = 931.49410242

#: neutron mass in MeV from CODATA 2018
m_n_MeV_C18 = 939.56542052

#: atomic mass unit in g from CODATA 2018
m_u_C18 = 1.6605390666e-27 * physical_constants.kilo

#: neutron mass in g from CODATA 2018
m_n_C18 = 1.67492749804e-27 * physical_constants.kilo
