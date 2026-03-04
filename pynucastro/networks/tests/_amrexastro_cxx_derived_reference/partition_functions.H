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


    // Fe52

    // this is log(partition function)

    inline AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 0, npts_1-1> Fe52_pf_array = {
        0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 3.999991999914312e-06, 2.1999758003515767e-05, 8.699621571943491e-05,
        0.00026096594542543384, 0.006963697081509575, 0.035510955889664755, 0.09259923578662213, 0.17220724583011637,
        0.26453344854718897, 0.3618510248890247, 0.45961295770683647, 0.5556743162453767, 0.7407490199034046,
        0.9195135330669838, 1.1002529420599954, 1.2949715212992456, 1.5184586008666647, 2.1077860146889784,
        2.9231615807191558, 3.893859034800475, 4.9344739331306915, 5.988961416889864, 7.038783541388542,
        8.06777619577889, 9.08704215563169, 10.094107912144779, 11.094375117313408, 13.559327152571864,
        15.990262279448435, 18.397412117013012, 20.790924485420227, 23.165612872315616, 25.527286881679668,
        27.88566333430213, 30.226275822885412, 32.55827480108575, 34.88236609930076, 37.1920183595179,
        39.494603452511946, 41.7901213782829, 44.07842051402947, 46.359186559628874, 48.63191868947324,
        50.895888946339504, 53.15830799760822, 55.41046223697537, 57.65993750465547, 59.90643313099847,
        62.14550481827019, 64.3824578963053, 66.61597196533687, 68.84562073247407, 71.07361272256216,
        73.2999973546708, 75.52327260920696, 77.74488863966732, 79.96498972270547, 84.40157534162978,
        88.83586272291235, 93.26655452956751, 97.69811509936366, 102.1291089050222, 106.56601751978464,
        111.00031864271057, 122.10466857715824,
    };

    constexpr amrex::Real Fe52_pf_threshold_T9 = 0.6;

    // Co55

    // this is log(partition function)

    inline AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 0, npts_1-1> Co55_pf_array = {
        0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 1.9999980000601777e-06, 3.299945551192844e-05, 0.0002149768908123338,
        0.0008776147834635239, 0.002650484360963482, 0.006510758935580232, 0.013783568660113304, 0.04542450669313041,
        0.11337689843044624, 0.2341190990271887, 0.4197354671252422, 0.6733637730790638, 1.3506671834767394,
        2.175887439948088, 3.068052935133617, 3.9982007016691985, 4.948759890378168, 5.910796644040527,
        6.887552571664617, 7.874739125171811, 8.869257522797287, 9.87302834505142, 12.392552212472792,
        14.927368073382404, 17.458346073576802, 19.980928412195695, 22.489707498190175, 24.983124837646084,
        27.463785196552635, 29.930601699902294, 32.384611307034916, 34.831446008873506, 37.26007182276292,
        39.68043995141335, 42.09363491595136, 44.49580471437484, 46.89326904555913, 49.27854585439829,
        51.66117365506587, 54.03786672021364, 56.40944468469464, 58.773587670688116, 61.13416002133251,
        63.491553350821555, 65.8460981827463, 68.19566348466455, 70.54312033183577, 72.89085982231036,
        75.23453177540692, 77.57658201060994, 79.91713370152783, 82.25818507534967, 86.9361146156202,
        91.61238072329202, 96.2900235580931, 100.96843290644959, 105.64972678791048, 110.33413387975575,
        115.02389413404445, 126.77320837707892,
    };

    constexpr amrex::Real Co55_pf_threshold_T9 = 1.5;

    // Ni56

    // this is log(partition function)

    inline AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 0, npts_1-1> Ni56_pf_array = {
        0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 9.999994999180668e-07, 1.7999838002017484e-05, 0.00014798904908051437,
        0.0006687763192560828, 0.0021007917908672197, 0.0051805577080849445, 0.0108093677159202, 0.033929817942853414,
        0.08125919821801637, 0.16557206393583876, 0.30308828495767837, 0.5109446166860524, 1.172482137234565,
        2.102913897864978, 3.1654750481410856, 4.272490747605575, 5.389071729816501, 6.498282149476434,
        7.605890001053122, 8.712759974960212, 9.814656338829513, 10.918718232265187, 13.664687668229691,
        16.38812278817138, 19.09362521719479, 21.774087461778173, 24.427033903554097, 27.056545465086103,
        29.664418719106976, 32.25599392921282, 34.82395533714435, 37.3778548584193, 39.91867374845115,
        42.44348030878443, 44.95807711272384, 47.46268883359117, 49.955692693611795, 52.438581179243556,
        54.91248475568723, 57.37829774665965, 59.83778360715438, 62.29201514356348, 64.73474686830077,
        67.18043280493549, 69.61404616033595, 72.04282585589065, 74.47118033617373, 76.89759077928012,
        79.31751257897871, 81.73470105471176, 84.14867938526342, 86.56428786666086, 91.38391256386228,
        96.2010760720766, 101.01398943805196, 105.82570189600392, 110.6374131490212, 115.4513381488714,
        120.26505308675254, 132.31206103765302,
    };

    constexpr amrex::Real Ni56_pf_threshold_T9 = 1.5;


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

    case Fe52:
        if (T9 > part_fun::Fe52_pf_threshold_T9) {
            part_fun::interpolate_pf(T9, pf_cache.index_temp_array_1, part_fun::temp_array_1, part_fun::Fe52_pf_array, logpf, dlogpf_dT9);
        }
        break;

    case Co55:
        if (T9 > part_fun::Co55_pf_threshold_T9) {
            part_fun::interpolate_pf(T9, pf_cache.index_temp_array_1, part_fun::temp_array_1, part_fun::Co55_pf_array, logpf, dlogpf_dT9);
        }
        break;

    case Ni56:
        if (T9 > part_fun::Ni56_pf_threshold_T9) {
            part_fun::interpolate_pf(T9, pf_cache.index_temp_array_1, part_fun::temp_array_1, part_fun::Ni56_pf_array, logpf, dlogpf_dT9);
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
        spec == He4 || spec == Fe52 || spec == Ni56
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
        spec == Co55
                       )
    {
        return 8.0_rt;
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
