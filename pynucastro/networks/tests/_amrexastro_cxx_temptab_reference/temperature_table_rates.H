#ifndef TEMPERATURE_TABLE_RATES_H
#define TEMPERATURE_TABLE_RATES_H

#include <AMReX_Array.H>

#include <tfactors.H>
#include <actual_network.H>


using namespace amrex::literals;

using namespace Rates;

// temperature / rate tabulation for O16 + p --> F17

namespace O16_p_to_F17_iliadis_data {
    inline AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 51> log_t9 = {
     -5.809142990314028  , -5.521460917862246  , -5.298317366548036  ,
     -5.115995809754082  , -4.961845129926823  , -4.8283137373023015 ,
     -4.710530701645918  , -4.605170185988091  , -4.509860006183766  ,
     -4.422848629194137  , -4.3428059215206005 , -4.268697949366879  ,
     -4.199705077879927  , -4.135166556742356  , -4.017383521085972  ,
     -3.912023005428146  , -3.6888794541139363 , -3.506557897319982  ,
     -3.2188758248682006 , -2.995732273553991  , -2.8134107167600364 ,
     -2.659260036932778  , -2.5257286443082556 , -2.4079456086518722 ,
     -2.3025850929940455 , -2.2072749131897207 , -2.120263536200091  ,
     -2.0402208285265546 , -1.9661128563728327 , -1.8971199848858813 ,
     -1.8325814637483102 , -1.7147984280919266 , -1.6094379124341003 ,
     -1.3862943611198906 , -1.2039728043259361 , -1.0498221244986778 ,
     -0.916290731874155  , -0.7985076962177716 , -0.6931471805599453 ,
     -0.5108256237659907 , -0.35667494393873245, -0.2231435513142097 ,
     -0.10536051565782628,  0.                 ,  0.22314355131420976,
     0.4054651081081644 ,  0.5596157879354227 ,  0.6931471805599453 ,
     0.9162907318741551 ,  1.0986122886681098 ,  1.252762968495368
    };

    inline AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 51> log_rate = {
     -92.93857540035118  , -82.5573056521023   , -75.17704808114306  ,
     -69.55076155001605  , -65.06026927373573  , -61.357534849389125 ,
     -58.22977030900149  , -55.537137188803406 , -53.184101816303524 ,
     -51.102378203514355 , -49.24174730308124  , -47.564466885054756 ,
     -46.04175152902774  , -44.648796232955824 , -42.18677046198662  ,
     -40.06826948447019  , -35.83012427872162  , -32.60450514439401  ,
     -27.91697135553557  , -24.595107052741724 , -22.06995494347394  ,
     -20.059032666551108 , -18.404807394796077 , -17.010181694135657 ,
     -15.811346515789314 , -14.766133968918934 , -13.842778985118468 ,
     -13.018453914013817 , -12.276709725129614 , -11.603835340330072 ,
     -10.989381068796463 ,  -9.905289174482757 ,  -8.974478048254198 ,
     -7.117735843770066 ,  -5.711109250042662 ,  -4.594230245949757 ,
     -3.6777417097034806,  -2.9069890637196023,  -2.246204759557938 ,
     -1.1628308609947609,  -0.3036759622061304,   0.4007875185570534,
     0.9928813340358675,   1.5000693753787353,   2.5087859235779435,
     3.2703291064571163,   3.872450230308323 ,   4.364753306003425 ,
     5.127529045567755 ,   5.694068782468562 ,   6.129268051255177
    };
}


namespace temp_tabular {

    template <int do_T_derivatives>
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    void rate_O16_p_to_F17_iliadis(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

        // O16 + p --> F17

        auto [_rate, _drate_dT] = interp_net::cubic_interp_uneven<do_T_derivatives>(
                                                   tfactors.lnT9,
                                                   O16_p_to_F17_iliadis_data::log_t9,
                                                   O16_p_to_F17_iliadis_data::log_rate);
        rate = std::exp(_rate);
        // we found dlog(rate)/dlog(T9)
        if constexpr (do_T_derivatives) {
            drate_dT = rate * tfactors.T9i * _drate_dT * 1.e-9;
        }
    }



    template <int do_T_derivatives, typename T>
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    void
    fill_rates(const tf_t& tfactors, T& rate_eval)
    {

        amrex::ignore_unused(tfactors);
        amrex::ignore_unused(rate_eval);

        [[maybe_unused]] amrex::Real rate;
        [[maybe_unused]] amrex::Real drate_dT;

        rate_O16_p_to_F17_iliadis<do_T_derivatives>(tfactors, rate, drate_dT);
        rate_eval.screened_rates(k_O16_p_to_F17_iliadis) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_O16_p_to_F17_iliadis) = drate_dT;

        }

    }


}
#endif
