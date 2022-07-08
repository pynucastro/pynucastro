#ifndef REACLIB_RATES_H
#define REACLIB_RATES_H

#include <AMReX.H>
#include <AMReX_Print.H>

#include <actual_network.H>

using namespace Rates;

struct rate_eval_t {
    Array1D<Real, 1, NumRates>  screened_rates;
    Array1D<Real, 1, NumRates>  dscreened_rates_dT;
    Array1D<Real, NrateReaclib+1, NrateReaclib+NrateTabular> add_energy_rate;
};

struct tf_t {
    Real T9;
    Real T9i;
    Real T943i;
    Real T923i;
    Real T913i;
    Real T913;
    Real T923;
    Real T953;
    Real lnT9;
};

AMREX_GPU_HOST_DEVICE AMREX_INLINE
tf_t evaluate_tfactors(const Real T)
{

    tf_t tf;
    tf.T9 = T / 1.e9_rt;
    tf.T9i = 1.0_rt / tf.T9;
    tf.T913 = std::cbrt(tf.T9);
    tf.T913i = 1.0_rt / tf.T913;
    tf.T923i = tf.T913i * tf.T913i;
    tf.T943i = tf.T9i * tf.T913i;
    tf.T923 = tf.T913 * tf.T913;
    tf.T953 = tf.T9 * tf.T923;
    tf.lnT9 = std::log(tf.T9);

    return tf;
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_c12_c12__he4_ne20(const tf_t& tfactors, Real& rate, Real& drate_dT) {

    // c12 + c12 --> he4 + ne20

    rate = 0.0;
    drate_dT = 0.0;

    Real ln_set_rate{0.0};
    Real dln_set_rate_dT9{0.0};
    Real set_rate{0.0};

    // cf88r
    ln_set_rate =  61.2863 + -84.165 * tfactors.T913i + -1.56627 * tfactors.T913
                         + -0.0736084 * tfactors.T9 + -0.072797 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    dln_set_rate_dT9 =  + -(1.0/3.0) * -84.165 * tfactors.T943i + (1.0/3.0) * -1.56627 * tfactors.T923i
                              + -0.0736084 + (5.0/3.0) * -0.072797 * tfactors.T923 + -0.666667 * tfactors.T9i;

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;

}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_c12_c12__n_mg23(const tf_t& tfactors, Real& rate, Real& drate_dT) {

    // c12 + c12 --> n + mg23

    rate = 0.0;
    drate_dT = 0.0;

    Real ln_set_rate{0.0};
    Real dln_set_rate_dT9{0.0};
    Real set_rate{0.0};

    // cf88r
    ln_set_rate =  -12.8056 + -30.1485 * tfactors.T9i + 11.4826 * tfactors.T913
                         + 1.82849 * tfactors.T9 + -0.34844 * tfactors.T953;

    dln_set_rate_dT9 =  30.1485 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 11.4826 * tfactors.T923i
                              + 1.82849 + (5.0/3.0) * -0.34844 * tfactors.T923;

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;

}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_c12_c12__p_na23(const tf_t& tfactors, Real& rate, Real& drate_dT) {

    // c12 + c12 --> p + na23

    rate = 0.0;
    drate_dT = 0.0;

    Real ln_set_rate{0.0};
    Real dln_set_rate_dT9{0.0};
    Real set_rate{0.0};

    // cf88r
    ln_set_rate =  60.9649 + -84.165 * tfactors.T913i + -1.4191 * tfactors.T913
                         + -0.114619 * tfactors.T9 + -0.070307 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    dln_set_rate_dT9 =  + -(1.0/3.0) * -84.165 * tfactors.T943i + (1.0/3.0) * -1.4191 * tfactors.T923i
                              + -0.114619 + (5.0/3.0) * -0.070307 * tfactors.T923 + -0.666667 * tfactors.T9i;

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;

}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_he4_c12__o16(const tf_t& tfactors, Real& rate, Real& drate_dT) {

    // c12 + he4 --> o16

    rate = 0.0;
    drate_dT = 0.0;

    Real ln_set_rate{0.0};
    Real dln_set_rate_dT9{0.0};
    Real set_rate{0.0};

    // nac2 
    ln_set_rate =  69.6526 + -1.39254 * tfactors.T9i + 58.9128 * tfactors.T913i + -148.273 * tfactors.T913
                         + 9.08324 * tfactors.T9 + -0.541041 * tfactors.T953 + 70.3554 * tfactors.lnT9;

    dln_set_rate_dT9 =  1.39254 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 58.9128 * tfactors.T943i + (1.0/3.0) * -148.273 * tfactors.T923i
                              + 9.08324 + (5.0/3.0) * -0.541041 * tfactors.T923 + 70.3554 * tfactors.T9i;

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;

    // nac2 
    ln_set_rate =  254.634 + -1.84097 * tfactors.T9i + 103.411 * tfactors.T913i + -420.567 * tfactors.T913
                         + 64.0874 * tfactors.T9 + -12.4624 * tfactors.T953 + 137.303 * tfactors.lnT9;

    dln_set_rate_dT9 =  1.84097 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 103.411 * tfactors.T943i + (1.0/3.0) * -420.567 * tfactors.T923i
                              + 64.0874 + (5.0/3.0) * -12.4624 * tfactors.T923 + 137.303 * tfactors.T9i;

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;

}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_n__p__weak__wc12(const tf_t& tfactors, Real& rate, Real& drate_dT) {

    // n --> p

    rate = 0.0;
    drate_dT = 0.0;

    Real ln_set_rate{0.0};
    Real dln_set_rate_dT9{0.0};
    Real set_rate{0.0};

    // wc12w
    ln_set_rate =  -6.78161;

    dln_set_rate_dT9 = 0.0;

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;

}



AMREX_GPU_HOST_DEVICE AMREX_INLINE
void
fill_reaclib_rates(const tf_t& tfactors, rate_eval_t& rate_eval)
{

    Real rate;
    Real drate_dT;

    rate_c12_c12__he4_ne20(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_c12_c12__he4_ne20) = rate;
    rate_eval.dscreened_rates_dT(k_c12_c12__he4_ne20) = drate_dT;

    rate_c12_c12__n_mg23(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_c12_c12__n_mg23) = rate;
    rate_eval.dscreened_rates_dT(k_c12_c12__n_mg23) = drate_dT;

    rate_c12_c12__p_na23(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_c12_c12__p_na23) = rate;
    rate_eval.dscreened_rates_dT(k_c12_c12__p_na23) = drate_dT;

    rate_he4_c12__o16(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_he4_c12__o16) = rate;
    rate_eval.dscreened_rates_dT(k_he4_c12__o16) = drate_dT;

    rate_n__p__weak__wc12(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_n__p__weak__wc12) = rate;
    rate_eval.dscreened_rates_dT(k_n__p__weak__wc12) = drate_dT;



}

#endif
