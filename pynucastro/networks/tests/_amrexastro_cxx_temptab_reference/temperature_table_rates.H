#ifndef TEMPERATURE_TABLE_RATES_H
#define TEMPERATURE_TABLE_RATES_H

#include <AMReX_Array.H>

#include <numbers>

#include <tfactors.H>
#include <actual_network.H>


using namespace amrex::literals;

using namespace Rates;

// temperature / rate tabulation for O16 + p --> F17

namespace O16_p_to_F17_iliadis_data {
    inline AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 51> log_t9 = {
     -2.5228787452803374 , -2.3979400086720375 , -2.3010299956639813 ,
     -2.221848749616356  , -2.154901959985743  , -2.0969100130080562 ,
     -2.0457574905606752 , -2.                 , -1.958607314841775  ,
     -1.9208187539523751 , -1.8860566476931633 , -1.853871964321762  ,
     -1.8239087409443189 , -1.7958800173440752 , -1.744727494896694  ,
     -1.6989700043360187 , -1.6020599913279623 , -1.5228787452803376 ,
     -1.3979400086720375 , -1.3010299956639813 , -1.2218487496163564 ,
     -1.154901959985743  , -1.0969100130080565 , -1.0457574905606752 ,
     -1.                 , -0.958607314841775  , -0.9208187539523752 ,
     -0.8860566476931632 , -0.8538719643217619 , -0.8239087409443188 ,
     -0.7958800173440752 , -0.744727494896694  , -0.6989700043360187 ,
     -0.6020599913279624 , -0.5228787452803376 , -0.4559319556497244 ,
     -0.3979400086720376 , -0.3467874862246563 , -0.3010299956639812 ,
     -0.22184874961635637, -0.1549019599857432 , -0.09691001300805639,
     -0.04575749056067512,  0.                 ,  0.09691001300805642,
     0.17609125905568124,  0.24303804868629442,  0.3010299956639812 ,
     0.3979400086720376 ,  0.47712125471966244,  0.5440680443502757
    };

    inline AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 51> log_rate = {
     -40.36271045232183   , -35.85418228550817   , -32.64897714741588   ,
     -30.20551195334083   , -28.255315936723115  , -26.647238808276168  ,
     -25.288867927693158  , -24.119472221801196  , -23.097561943801335  ,
     -22.193480865919295  , -21.385419133002515  , -20.65698550284923   ,
     -19.99567862621736   , -19.390725827595414  , -18.321481620959887  ,
     -17.40142833651786   , -15.560825260156532  , -14.159956669396506  ,
     -12.124186611160242  , -10.681519274825483  ,  -9.584859647804127  ,
     -8.711527199400217  ,  -7.9931062920521    ,  -7.387428045934824  ,
     -6.866780543267506  ,  -6.412850501745656  ,  -6.011842527443247  ,
     -5.653842697767992  ,  -5.331707289551779  ,  -5.039481657219292  ,
     -4.772627557710364  ,  -4.301812430133878  ,  -3.897566294318664  ,
     -3.091193400594326  ,  -2.480303232840147  ,  -1.9952488444089989 ,
     -1.5972229303896526 ,  -1.262489309326524  ,  -0.9755143323008331 ,
     -0.5050110263168319 ,  -0.13188479467278288,   0.17405980772502544,
     0.4312028845565166 ,   0.6514718521990425 ,   1.089551882886454  ,
     1.420285884941918  ,   1.6817837664678814 ,   1.8955882756662628 ,
     2.2268575702887237 ,   2.472902651803664  ,   2.661907292766021
    };
}


namespace temp_tabular {

    constexpr amrex::Real ln10_inv = 1.0_rt / std::numbers::ln10_v<amrex::Real>;

    template <int do_T_derivatives>
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    void rate_O16_p_to_F17_iliadis(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

        // O16 + p --> F17

        amrex::Real log_t9 = tfactors.lnT9 * ln10_inv;
        auto [_rate, _drate_dT] = interp_net::cubic_interp_uneven<do_T_derivatives>(
                                                   log_t9,
                                                   O16_p_to_F17_iliadis_data::log_t9,
                                                   O16_p_to_F17_iliadis_data::log_rate);
        rate = amrex::Math::exp10(_rate);
        // we found dlog10(rate)/dlog10(T9)
        if constexpr (do_T_derivatives) {
            drate_dT = (rate / tfactors.T9) * _drate_dT * 1.e-9;
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
