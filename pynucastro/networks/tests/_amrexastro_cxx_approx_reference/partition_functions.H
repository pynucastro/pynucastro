#ifndef PARTITION_FUNCTIONS_H
#define PARTITION_FUNCTIONS_H

#include <AMReX_REAL.H>
#include <AMReX_Array.H>
#include <AMReX_Math.H>

#include <tfactors.H>
#include <fundamental_constants.H>
#include <network_properties.H>
#include <interp_tools.H>

using namespace amrex::literals;
using namespace Species;

namespace part_fun {


    // interpolation routine

    template <typename T>
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    void interpolate_pf(const amrex::Real t9, int idx,
                        const T& temp_array, const T& pf_array,
                        amrex::Real& pf, amrex::Real& dpf_dT) {

        // if the index is -1, we are either outside of the bounds or
        // we never computed it.  try recomputing here
        if (idx == -1) {
            idx = interp_net::find_index(t9, temp_array);
        }

        // if the index is still -1, then we are out-of-bounds, so
        // just skip the interpolation

        if (idx != -1) {

            // we have idx such that:
            // temp_array[idx] <= t9 < temp_array[idx+1]
            AMREX_ASSERT(t9 >= temp_array(idx) && t9 < temp_array(idx+1));

            // construct the slope -- this is
            // (log10(pf_{i+1}) - log10(pf_i)) / (T_{i+1} - T_i)

            amrex::Real slope = (pf_array(idx+1) - pf_array(idx)) /
                                (temp_array(idx+1) - temp_array(idx));

            // find the PF

            amrex::Real log10_pf = pf_array(idx) + slope * (t9 - temp_array(idx));
            pf = amrex::Math::exp10(log10_pf);

            // find the derivative (with respect to T, not T9)

            amrex::Real dpf_dT9 = pf * M_LN10 * slope;
            dpf_dT = dpf_dT9 * 1.e-9_rt;

        } else {

            // T < the smallest T or >= the largest T in the partition function table
            pf = 1.0;
            dpf_dT = 0.0;

        }

    }

    struct pf_cache_t {
        // Store the coefficient and derivative adjacent in memory, as they're
        // always accessed at the same time.
        // The entries will be default-initialized to zero, which is fine since
        // log10(x) is never zero.
        amrex::Array2D<amrex::Real, 1, NumSpecTotal, 1, 2, amrex::Order::C> data{};
    };

}

// get the partition function for nucleus inuc.  Here pf_cache
// should already contain the index into the temp / pf arrays

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void get_partition_function(const int inuc,
                            [[maybe_unused]] const tf_t& tfactors,
                            part_fun::pf_cache_t& pf_cache,
                            amrex::Real& pf, amrex::Real& dpf_dT) {

    pf = 1.0_rt;
    dpf_dT = 0.0_rt;

    // inuc is the 1-based index for the species

    switch (inuc) {


    default:

        pf = 1.0_rt;
        dpf_dT = 0.0_rt;

    }

}


// a version of get_parition_function where we have not cached the
// index into the arrays

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void get_partition_function(const int inuc,
                            [[maybe_unused]] const tf_t& tfactors,
                            amrex::Real& pf, amrex::Real& dpf_dT) {

    // create a dummy cache.  This will initialize the indices
    // to -1 and we will force a search for the correct index
    // when we do the interpolation.

    part_fun::pf_cache_t pf_cache;

    get_partition_function(inuc, tfactors, pf_cache, pf, dpf_dT);
}


AMREX_GPU_HOST_DEVICE AMREX_INLINE
void get_partition_function_cached(const int inuc, const tf_t& tfactors,
                                   part_fun::pf_cache_t& pf_cache,
                                   amrex::Real& pf, amrex::Real& dpf_dT) {
    if (pf_cache.data(inuc, 1) != 0.0_rt) {
        // present in cache
        amrex::ignore_unused(tfactors);
        pf = pf_cache.data(inuc, 1);
        dpf_dT = pf_cache.data(inuc, 2);
    } else {
        get_partition_function(inuc, tfactors, pf_cache, pf, dpf_dT);
        pf_cache.data(inuc, 1) = pf;
        pf_cache.data(inuc, 2) = dpf_dT;
    }
}

// spins

AMREX_GPU_HOST_DEVICE AMREX_INLINE
constexpr amrex::Real get_spin_state(const int inuc) {

    amrex::Real spin = -1.0;

    switch (inuc) {  // NOLINT(bugprone-switch-missing-default-case)

    case He4:
    case Mg24:
    case Si28:
    case S32:
        spin = 1;
        break;

    case P31:
        spin = 2;
        break;

    case Al27:
        spin = 6;
        break;


    }

    return spin;

}


#endif
