#include <actual_network.H>
#include <amrex_bridge.H>

namespace network
{
    Array1D<Real, 1, NumSpec> bion;
    Array1D<Real, 1, NumSpec> mion;
}

void actual_network_init()
{
    using namespace Species;
    using namespace network;

    // binding energies per nucleon in MeV
    Array1D<Real, 1, NumSpec> ebind_per_nucleon;

    ebind_per_nucleon(N) = 0.0_rt;
    ebind_per_nucleon(H1) = 0.0_rt;
    ebind_per_nucleon(He4) = 7.073915614499924_rt;
    ebind_per_nucleon(C12) = 7.680144581999836_rt;
    ebind_per_nucleon(O16) = 7.97620721324995_rt;
    ebind_per_nucleon(Ne20) = 8.032241192000038_rt;
    ebind_per_nucleon(Na23) = 8.111493582782714_rt;
    ebind_per_nucleon(Mg23) = 7.90112268991304_rt;

    // convert to binding energies per nucleus in MeV
    for (int i = 1; i <= NumSpec; ++i) {
        bion(i) = ebind_per_nucleon(i) * aion[i-1];
    }

    // Set the mass -- this will be in grams
    mion(N) = 1.674927498034172e-24_rt;
    mion(H1) = 1.6735328377636005e-24_rt;
    mion(He4) = 6.646479071584587e-24_rt;
    mion(C12) = 1.99264687992e-23_rt;
    mion(O16) = 2.6560180592333686e-23_rt;
    mion(Ne20) = 3.3198227947612416e-23_rt;
    mion(Na23) = 3.817541002484691e-23_rt;
    mion(Mg23) = 3.8182640828719474e-23_rt;

}
