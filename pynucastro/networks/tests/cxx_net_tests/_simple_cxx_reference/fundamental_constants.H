#ifndef FUNDAMENTAL_CONSTANTS_H
#define FUNDAMENTAL_CONSTANTS_H

// Constants are obtained from scipy.constants.

#include <amrex_bridge.H>

namespace C
{
    // speed of light in vacuum
    constexpr Real c_light = 29979245800.0;  // cm/s

    // boltzmann's constant
    constexpr Real k_B = 1.380649e-16;  // erg/K

    // planck's constant over 2pi
    constexpr Real hbar = 1.0545718176461563e-27;  // erg

    // avogradro's Number
    constexpr Real n_A = 6.02214076e+23;  // mol^-1

    // convert eV to erg
    constexpr Real ev2erg = 1.602176634e-12;

    // convert MeV to eV
    constexpr Real MeV2eV = 1.0e6;

    // convert MeV to grams
    constexpr Real MeV2gr  = (MeV2eV * ev2erg) / (c_light * c_light);

    // conversion factor for nuclear energy generation rate
    constexpr Real enuc_conv2 = -n_A * c_light * c_light;

    // mass of proton
    constexpr Real m_p = 1.67262192595e-24;  // g

    // mass of neutron
    constexpr Real m_n = 1.6749275005600003e-24;  // g

    // mass of electron
    constexpr Real m_e = 9.1093837139e-28;  // g

    // atomic mass unit
    constexpr Real m_u = 1.66053906892e-24; // g

    // electron charge
    // NIST: q_e = 1.602176565e-19 C
    //
    // C is the SI unit Coulomb; in cgs we have the definition:
    //     1 C = 0.1 * |c_light| * 1 statC
    // where statC is the cgs unit statCoulomb; 1 statC = 1 erg^1/2 cm^1/2
    // and |c_light| is the speed of light in cgs (but without units)
    constexpr Real q_e = 4.803204712570263e-10;  // erg^1/2 cm^1/2

}

#endif
