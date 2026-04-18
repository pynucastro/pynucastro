#ifndef FUNDAMENTAL_CONSTANTS_H
#define FUNDAMENTAL_CONSTANTS_H

// Fundamental constants taken from NIST's 2010 CODATA recommended values

#include <amrex_bridge.H>

namespace C
{
    // speed of light in vacuum
    constexpr Real c_light = 2.99792458e10;  // cm/s

    // boltzmann's constant
    constexpr Real k_B = 1.3806488e-16;  // erg/K

    // planck's constant over 2pi
    constexpr Real hbar = 1.054571726e-27;  // erg s

    // avogradro's Number
    constexpr Real n_A = 6.02214129e23;  // mol^-1

    // convert eV to erg
    constexpr Real ev2erg = 1.602176487e-12;

    // convert MeV to eV
    constexpr Real MeV2eV = 1.0e6;

    // convert MeV to grams
    constexpr Real MeV2gr  = (MeV2eV * ev2erg) / (c_light * c_light);

    // conversion factor for nuclear energy generation rate
    constexpr Real enuc_conv2 = -n_A * c_light * c_light;

    // mass of proton
    constexpr Real m_p = 1.672621777e-24;  // g

    // mass of neutron
    constexpr Real m_n = 1.674927351e-24;  // g

    // mass of electron
    constexpr Real m_e = 9.10938291e-28;  // g

    // atomic mass unit
    constexpr Real m_u = 1.6605390666e-24; // g

    // electron charge
    // NIST: q_e = 1.602176565e-19 C
    //
    // C is the SI unit Coulomb; in cgs we have the definition:
    //     1 C = 0.1 * |c_light| * 1 statC
    // where statC is the cgs unit statCoulomb; 1 statC = 1 erg^1/2 cm^1/2
    // and |c_light| is the speed of light in cgs (but without units)
    constexpr Real q_e = 4.80320451e-10;  // erg^1/2 cm^1/2

}

#endif
