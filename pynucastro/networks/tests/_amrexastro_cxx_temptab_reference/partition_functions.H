#ifndef PARTITION_FUNCTIONS_H
#define PARTITION_FUNCTIONS_H

#include <numbers>

#include <AMReX_REAL.H>
#include <AMReX_Array.H>
#include <AMReX_Math.H>
#include <AMReX_ConstexprFor.H>

#include <tfactors.H>
#include <fundamental_constants.H>
#include <network_properties.H>
#include <interp_tools.H>

using namespace amrex::literals;
using namespace Species;

namespace part_fun {

    constexpr int npts_1 = 72;

    // this is T9

    inline AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 0, npts_1-1> temp_array_1 = {
        0.01, 0.15, 0.2, 0.3, 0.4,
        0.5, 0.6, 0.7, 0.8, 0.9,
        1.0, 1.5, 2.0, 2.5, 3.0,
        3.5, 4.0, 4.5, 5.0, 6.0,
        7.0, 8.0, 9.0, 10.0, 12.0,
        14.0, 16.0, 18.0, 20.0, 22.0,
        24.0, 26.0, 28.0, 30.0, 35.0,
        40.0, 45.0, 50.0, 55.0, 60.0,
        65.0, 70.0, 75.0, 80.0, 85.0,
        90.0, 95.0, 100.0, 105.0, 110.0,
        115.0, 120.0, 125.0, 130.0, 135.0,
        140.0, 145.0, 150.0, 155.0, 160.0,
        165.0, 170.0, 175.0, 180.0, 190.0,
        200.0, 210.0, 220.0, 230.0, 240.0,
        250.0, 275.0,
    };


    // O16

    // this is log(partition function)

    inline AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 0, npts_1-1> O16_pf_array = {
        0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.02955880224154443,
        0.08617769624105241, 0.16551443847757333, 0.27763173659827955, 0.4252677354043441, 0.6043159668533296,
        0.7975071958841882, 1.007957920399979, 1.2296405510745139, 1.4609379041156563, 2.0681277817795625,
        2.70805020110221, 3.3775875160230218, 4.062165663857866, 4.7535901911063645, 5.4510384535657,
        6.142037405587356, 6.828712071641684, 7.512071245835466, 8.185907481482324, 8.85366542803745,
        9.517825071724143, 10.177324218165856, 10.831706855275558, 11.48246625748552, 12.128111104060462,
        12.774223335915433, 13.41652441595382, 14.054527458434775, 14.690979295318174, 15.324022551808413,
        15.956752500549557, 16.588099280204055, 17.216707939626428, 17.842646370492925, 18.469470908121796,
        19.09362521719479, 19.715407911546766, 20.337603356134426, 20.954377557909798, 22.19802884605391,
        23.43796058076729, 24.670655986211848, 25.905049387238495, 27.136724794113768, 28.368185181905268,
        29.602320498988682, 32.68726692127686,
    };

    constexpr amrex::Real O16_pf_threshold_T9 = 10.0;

    // F17

    // this is log(partition function)

    inline AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 0, npts_1-1> F17_pf_array = {
        0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.009950330853168092, 0.01980262729617973, 0.02955880224154443, 0.04879016416943205,
        0.058268908123975824, 0.0769610411361284, 0.08617769624105241, 0.09531017980432493, 0.12221763272424911,
        0.13976194237515863, 0.1570037488096647, 0.17395330712343798, 0.19885085874516517, 0.26236426446749106,
        0.3293037471426003, 0.41871033485818504, 0.5247285289349821, 0.6471032420585384, 0.7884573603642703,
        0.9400072584914712, 1.1019400787607843, 1.2781522025001875, 1.4678743481123135, 1.9810014688665833,
        2.5416019934645457, 3.1441522786722644, 3.7612001156935624, 4.388257184424518, 5.0238805208462765,
        5.655991810819852, 6.2878585601617845, 6.917705609835305, 7.549609165154532, 8.183118079394745,
        8.814330422638774, 9.441452092939569, 10.077440859659566, 10.707728780601661, 11.33857207782545,
        11.970350312009105, 12.601487417784837, 13.233904752137237, 13.864300722133706, 14.498607402670718,
        15.131918791619999, 15.765697263786848, 16.3957273875566, 17.034386382832476, 17.667783558986645,
        18.301897207962398, 18.939474537367534, 19.571252771551187, 20.207427671356875, 21.479387816667746,
        22.755353682242777, 24.033808850340435, 25.311289864099532, 26.592562750080187, 27.877881193860073,
        29.165735482166713, 32.393195050726305,
    };

    constexpr amrex::Real F17_pf_threshold_T9 = 1.0;


    // interpolation routine

    template <typename T>
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    void interpolate_pf(const amrex::Real t9, int idx,
                        const T& temp_array, const T& pf_array,
                        amrex::Real& logpf, amrex::Real& dlogpf_dT9) {

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

            // find the derivative: dlog(pf)/dT9 -- this is
            // (log(pf_{i+1}) - log(pf_i)) / (T_{i+1} - T_i)

            dlogpf_dT9 = (pf_array(idx+1) - pf_array(idx)) /
                           (temp_array(idx+1) - temp_array(idx));

            // find log(PF)

            logpf = pf_array(idx) + dlogpf_dT9 * (t9 - temp_array(idx));

        } else {

            // T < the smallest T or >= the largest T in the partition function table
            logpf = 0.0;
            dlogpf_dT9 = 0.0;

        }

    }

    struct pf_cache_t {
        // Store the coefficient and derivative adjacent in memory, as they're
        // always accessed at the same time.
        amrex::Array2D<amrex::Real, 1, NumSpecTotal, 1, 2, amrex::Order::C> data;
        int index_temp_array_1{-1};

        AMREX_GPU_HOST_DEVICE AMREX_INLINE
        pf_cache_t() {
            // The entries will be default-initialized to -900,
            // log(partition_function) should never be -900.
            for (int i = 1; i <= NumSpecTotal; ++i) {
                    data(i, 1) = -900.0_rt;
                    data(i, 2) = -900.0_rt;
            }
        }
    };

}

// get the partition function for nucleus inuc.  Here pf_cache
// should already contain the index into the temp / pf arrays

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void get_partition_function(const int inuc,
                            [[maybe_unused]] const amrex::Real T9,
                            const part_fun::pf_cache_t& pf_cache,
                            amrex::Real& logpf, amrex::Real& dlogpf_dT9) {

    logpf = 0.0_rt;
    dlogpf_dT9 = 0.0_rt;

    // inuc is the 1-based index for the species

    switch (inuc) {

    case O16:
        if (T9 > part_fun::O16_pf_threshold_T9) {
            part_fun::interpolate_pf(T9, pf_cache.index_temp_array_1, part_fun::temp_array_1, part_fun::O16_pf_array, logpf, dlogpf_dT9);
        }
        break;

    case F17:
        if (T9 > part_fun::F17_pf_threshold_T9) {
            part_fun::interpolate_pf(T9, pf_cache.index_temp_array_1, part_fun::temp_array_1, part_fun::F17_pf_array, logpf, dlogpf_dT9);
        }
        break;


    default:

        logpf = 0.0_rt;
        dlogpf_dT9 = 0.0_rt;

    }

}


// a version of get_partition_function where we have not cached the
// index into the arrays. Uses tfactors instead of T9 as argument,
// and returns pf and dpf_dT instead of logpf and dlogpf_dT9

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void get_partition_function(const int inuc,
                            [[maybe_unused]] const tf_t& tfactors,
                            amrex::Real& pf, amrex::Real& dpf_dT) {

    // create a dummy cache.  This will initialize the indices
    // to -1 and we will force a search for the correct index
    // when we do the interpolation.

    part_fun::pf_cache_t pf_cache;

    amrex::Real logpf, dlogpf_dT9;
    get_partition_function(inuc, tfactors.T9, pf_cache, logpf, dlogpf_dT9);

    // Convert logpf to pf and dlogpf_dT9 to dpf_dT
    pf = std::exp(logpf);
    dpf_dT = pf * dlogpf_dT9 * 1.e-9_rt;
}


AMREX_GPU_HOST_DEVICE AMREX_INLINE
void get_partition_function_cached(const int inuc, const amrex::Real T9,
                                   part_fun::pf_cache_t& pf_cache,
                                   amrex::Real& logpf, amrex::Real& dlogpf_dT9) {
    if (pf_cache.data(inuc, 1) != -900.0_rt) {
        // present in cache
        amrex::ignore_unused(T9);
        logpf = pf_cache.data(inuc, 1);
        dlogpf_dT9 = pf_cache.data(inuc, 2);
    } else {
        get_partition_function(inuc, T9, pf_cache, logpf, dlogpf_dT9);
        pf_cache.data(inuc, 1) = logpf;
        pf_cache.data(inuc, 2) = dlogpf_dT9;
    }
}


template <int spec>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
constexpr amrex::Real get_spin_state() {

    static_assert(spec >= 1 && spec <= NumSpec);

    // Set the spin states (2J + 1) of the element

    if constexpr (
        spec == O16
                  )
    {
        return 1.0_rt;
    }
    else if constexpr (
        spec == H1
                       )
    {
        return 2.0_rt;
    }
    else if constexpr (
        spec == F17
                       )
    {
        return 6.0_rt;
    }

    // Return -1 if we don't recognize the species.
    return -1.0_rt;
}


// Legacy (non-templated) interface

AMREX_GPU_HOST_DEVICE AMREX_INLINE
amrex::Real get_spin_state (int spec) {

    amrex::Real s = -1.0_rt;

    // Set the spin states (2J + 1) of the element

    amrex::constexpr_for<1, NumSpec+1>([&] (auto n) {
        if (n == spec) {
            s = get_spin_state<n>();
        }
    });

    return s;
}


#endif
