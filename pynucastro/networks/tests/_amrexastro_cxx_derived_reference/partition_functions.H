#ifndef PARTITION_FUNCTIONS_H
#define PARTITION_FUNCTIONS_H

#include <AMReX_REAL.H>
#include <AMReX_Array.H>

#include <tfactors.H>
#include <fundamental_constants.H>
#include <network_properties.H>

using namespace amrex::literals;
using namespace Species;

namespace part_fun {

    constexpr int npts_1 = 72;

    // this is T9

    extern AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 0, npts_1> temp_array_1;

    // this is log10(partition function)

    extern AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 0, npts_1> Fe52_pf_array;

    // this is log10(partition function)

    extern AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 0, npts_1> Co55_pf_array;

    // this is log10(partition function)

    extern AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 0, npts_1> Ni56_pf_array;



    // interpolation routine

    template <typename T>
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    void interpolate_pf(const amrex::Real t9, const T& temp_array, const T& pf_array,
                        amrex::Real& pf, amrex::Real& dpf_dT) {

        if (t9 >= temp_array.lo() && t9 < temp_array.hi()) {

            // find the largest temperature element <= t9 using a binary search

            int left = temp_array.lo();
            int right = temp_array.hi();

            while (left < right) {
                int mid = (left + right) / 2;
                if (temp_array(mid) > t9) {
                    right = mid;
                } else {
                    left = mid + 1;
                }
            }

            const int idx = right - 1;

            // now we have temp_array[idx] <= t9 < temp_array[idx+1]

            // construct the slope -- this is (log10(pf_{i+1}) - log10(pf_i)) / (T_{i+1} - T_i)

            amrex::Real slope = (pf_array(idx+1) - pf_array(idx)) /
                                (temp_array(idx+1) - temp_array(idx));

            // find the PF

            amrex::Real log10_pf = pf_array(idx) + slope * (t9 - temp_array(idx));
            pf = std::pow(10.0_rt, log10_pf);

            // find the derivative (with respect to T, not T9)

            amrex::Real dpf_dT9 = pf * M_LN10 * slope;
            dpf_dT = dpf_dT9 / 1.e9_rt;

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

// main interface

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void get_partition_function(const int inuc, [[maybe_unused]] const tf_t& tfactors,
                            amrex::Real& pf, amrex::Real& dpf_dT) {

    // inuc is the 1-based index for the species

    switch (inuc) {

    case Fe52:
        part_fun::interpolate_pf(tfactors.T9, part_fun::temp_array_1, part_fun::Fe52_pf_array, pf, dpf_dT);
        break;

    case Co55:
        part_fun::interpolate_pf(tfactors.T9, part_fun::temp_array_1, part_fun::Co55_pf_array, pf, dpf_dT);
        break;

    case Ni56:
        part_fun::interpolate_pf(tfactors.T9, part_fun::temp_array_1, part_fun::Ni56_pf_array, pf, dpf_dT);
        break;


    default:

        pf = 1.0_rt;
        dpf_dT = 0.0_rt;

    }

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
        get_partition_function(inuc, tfactors, pf, dpf_dT);
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
    case Fe52:
    case Ni56:
        spin = 1;
        break;

    case H1:
        spin = 2;
        break;

    case Co55:
        spin = 8;
        break;


    }

    return spin;

}


#endif
