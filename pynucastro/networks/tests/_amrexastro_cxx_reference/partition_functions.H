#ifndef PARTITION_FUNCTIONS_H
#define PARTITION_FUNCTIONS_H

#include <AMReX_REAL.H>
#include <AMReX_Array.H>

#include <tfactors.H>
#include <fundamental_constants.H>
#include <network_properties.H>

using namespace amrex;
using namespace Species;

namespace part_fun {



    // interpolation routine

    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    void interpolate_pf(const Real t9, const int npts, const Real* temp_array, const Real* pf_array,
                        Real& pf, Real& dpf_dT) {

        // find the index of the first temperature element we are larger than

        int idx = -1;

        for (int i = 0; i < npts-1; ++i) {
            if (t9 >= temp_array[i] && t9 < temp_array[i+1]) {
                idx = i;
                break;
            }
        }

        if (idx >= 0) {

            // construct the slope -- this is (log10(pf_{i+1}) - log10(pf_i)) / (T_{i+1} - T_i)

            Real slope = (pf_array[idx+1] - pf_array[idx]) / (temp_array[idx+1] - temp_array[idx]);

            // find the PF

            Real log10_pf = pf_array[idx] + slope * (t9 - temp_array[idx]);
            pf = std::pow(10.0_rt, log10_pf);

            // find the derivative (with respect to T, not T9)

            Real dpf_dT9 = pf * std::log(10.0_rt) * slope;
            dpf_dT = dpf_dT9 / 1.e9_rt;

        } else {

            // T < the smallest T in the partition function table
            pf = 1.0;
            dpf_dT = 0.0;

        }

    }

}

// main interface

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void get_partition_function(const int inuc, [[maybe_unused]] const tf_t& tfactors,
                            Real& pf, Real& dpf_dT) {

    // inuc is the 1-based index for the species

    switch (inuc) {


    default:

        pf = 1.0_rt;
        dpf_dT = 0.0_rt;

    }

}

// spins

AMREX_GPU_HOST_DEVICE AMREX_INLINE
constexpr Real get_spin_state(const int inuc) {

    Real spin = -1.0;

    switch (inuc) {

    case N:
        spin = 2;
        break;

    case H1:
        spin = 2;
        break;

    case He4:
        spin = 1;
        break;

    case C12:
        spin = 1;
        break;

    case O16:
        spin = 1;
        break;

    case Ne20:
        spin = 1;
        break;

    case Ne23:
        spin = 6;
        break;

    case Na23:
        spin = 4;
        break;

    case Mg23:
        spin = 4;
        break;


    }

    return spin;

}


#endif
