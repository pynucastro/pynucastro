#ifndef SCREEN_H
#define SCREEN_H

#include <cmath>
#include <iostream>
#include <numbers>
#include <string>

#include <amrex_bridge.H>

#include <network_properties.H>
#include <fundamental_constants.H>
#include <screen_data.H>
#include <autodiff.H>

// Composition partials at fixed temperature, density, and reacting pair.
// Ye = sum(Z_i Y_i), S2 = sum(Z_i^2 Y_i).
using screening_dual_t = dual::Dual<Real, dual::DerivativePair<Real>>;

template <typename number_t = Real>
struct plasma_state_t {
    number_t temp = 0.0_rt;
    number_t zbar = 0.0_rt;
    number_t z2bar = 0.0_rt;
    number_t n_e = 0.0_rt;
    number_t gamma_e_fac = 0.0_rt;
};


template <typename number_t>
inline
void
fill_plasma_state(plasma_state_t<number_t>& state, const Real& temp,
                  const Real dens, Array1D<Real, 1, NumSpec> const& y) {

    using std::cbrt;

    Real sum = 0.0_rt;
    for (int n = 1; n <= NumSpec; n++) {
        sum += y(n);
    }
    Real abar = 1.0_rt / sum;
    Real ytot = sum;

    number_t ye = 0.0_rt;
    number_t s2 = 0.0_rt;
    for (int n = 1; n <= NumSpec; n++) {
        ye = ye + zion[n-1]*y(n);
        s2 = s2 + zion[n-1]*zion[n-1]*y(n);
    }
    if constexpr (std::is_same_v<number_t, screening_dual_t>) {
        // abar cancels from the screening expressions, so hold it fixed.
        ye.seed({1.0_rt, 0.0_rt});
        s2.seed({0.0_rt, 1.0_rt});
    }

    state.zbar = ye * abar;
    state.z2bar = s2 * abar;

    // ntot
    Real rr = dens * ytot;

    state.temp = temp;

    // Electron number density
    // zbar * ntot works out to sum(z[i] * n[i]), after cancelling terms
    state.n_e = state.zbar * rr * C::n_A;

    // precomputed part of Gamma_e, from Chugunov 2009 eq. 6
    Real gamma_e_constants =
        C::q_e*C::q_e/C::k_B * std::pow(4.0_rt/3.0_rt*std::numbers::pi, 1.0_rt/3.0_rt);
    state.gamma_e_fac = gamma_e_constants * cbrt(state.n_e);
}


template <typename number_t>
inline
number_t chugunov2007 (const plasma_state_t<number_t>& state,
                   const scrn::screen_factors_t& scn_fac)
{
    using std::cbrt;
    using std::cos;
    using std::max;
    using std::pow;
    using std::sqrt;

    // Calculates screening factors based on Chugunov et al. 2007, following the
    // the approach in Yakovlev 2006 to extend to a multi-component plasma.
    //
    // References:
    //   Chugunov, DeWitt, and Yakovlev 2007, PhRvD, 76, 025028
    //   Yakovlev, Gasques, Afanasjev, Beard, and Wiescher 2006, PhRvC, 74, 035803
    //   Chugunov and DeWitt 2009, PhRvC, 80, 014611

    // input:
    // state   = plasma state (T, rho, abar, zbar, etc.)
    // scn_fac = screening factors for A and Z

    // Plasma temperature T_p
    // This formula comes from working backwards from zeta_ij (Chugunov 2009 eq. 12)
    // through Chugunov 2007 eq. 3 to Chugunov 2007 eq. 2.
    // Ultimately, the changes from the expression in Chugunov 2007 are:
    //   Z^2 -> Z1 * Z2
    //   n_i -> n_e / ztilde^3, where ztilde = (Z1^(1/3) + Z2^(1/3)) / 2
    //   m_i -> 2 mu12 (reduced mass)
    // This prescription reduces to the expressions from Chugunov 2007 in the case
    // of an OCP, and to Chugunov 2009 in the case of a binary ionic mixture.
    // This also matches Yakovlev et al. 2006, eq. 10.
    //
    // For reference, MESA r21.12.1 does:
    //   Z^2 -> Z1 * Z2
    //   n_i -> n_e / zbar (=ntot)
    //   m_i -> m_u * abar
    // Sam Jones' Fortran implementation (linked in issue #146) does:
    //   Z^2 -> zbar^2
    //   n_i -> ntot
    //   m_i -> m_u * abar
    Real mu12 = scn_fac.a1 * scn_fac.a2 / (scn_fac.a1 + scn_fac.a2);
    Real z_factor = scn_fac.z1 * scn_fac.z2;
    number_t n_i = state.n_e / scn_fac.ztilde3;
    Real m_i = 2.0_rt * mu12 / C::n_A;

    const Real T_p_factor = C::hbar/C::k_B*C::q_e * sqrt(4.0_rt*std::numbers::pi);
    number_t T_p = T_p_factor * sqrt(z_factor * n_i / m_i);

    // Normalized temperature
    number_t T_norm = state.temp / T_p;

    // The fit has only been verified down to T ~ 0.1 T_p, below which the rate
    // should be nearly temperature-independent (in the pycnonuclear regime),
    // and we clip the temperature to 0.1 T_p at small T.
    // start the transition here
    constexpr Real T_norm_fade = 0.2_rt;
    constexpr Real T_norm_min = 0.1_rt;

    if (T_norm < T_norm_min) {
        // clip temperature to the minimum value
        T_norm = T_norm_min;
    } else if (T_norm <= T_norm_fade) {
        // blend using a cosine, after MESA
        constexpr Real delta_T = T_norm_fade - T_norm_min;
        number_t tmp = std::numbers::pi * (T_norm - T_norm_min) / delta_T;
        number_t f = 0.5_rt * (1.0_rt - cos(tmp));
        T_norm = (1.0_rt - f) * T_norm_min + f * T_norm;
    }
    number_t inv_T_norm = 1.0_rt / T_norm;

    // Coulomb coupling parameter from Yakovlev 2006 eq. 10
    number_t Gamma = state.gamma_e_fac*scn_fac.z1*scn_fac.z2 / (scn_fac.ztilde*T_p) * inv_T_norm;

    // The fit for Gamma is only applicable up to ~600, so smoothly cap its value
    constexpr Real Gamma_fade = 590;
    constexpr Real Gamma_max = 600;

    if (Gamma > Gamma_max) {
        // clip Gamma to the max value
        Gamma = Gamma_max;
    } else if (Gamma >= Gamma_fade) {
        // blend using a cosine, after MESA
        constexpr Real delta_gamma = Gamma_max - Gamma_fade;
        number_t tmp = std::numbers::pi * (Gamma - Gamma_fade) / delta_gamma;
        number_t f = 0.5_rt * (1.0_rt - cos(tmp));
        Gamma = (1.0_rt - f) * Gamma + f * Gamma_max;
    }

    // Chugunov 2007 eq. 3
    constexpr Real zeta_factor = 4.0_rt / (3.0_rt * std::numbers::pi * std::numbers::pi);
    number_t zeta = cbrt(zeta_factor * (inv_T_norm*inv_T_norm));

    // Gamma tilde from Chugunov 2007 eq. 21
    constexpr Real fit_alpha = 0.022_rt;
    number_t fit_beta = 0.41_rt - 0.6_rt / Gamma;
    number_t fit_gamma = 0.06_rt + 2.2_rt / Gamma;
    // Polynomial term in Gamma tilde
    number_t poly = 1.0_rt + zeta*(fit_alpha + zeta*(fit_beta + fit_gamma*zeta));

    number_t gamtilde = Gamma / cbrt(poly);

    // fit parameters just after Chugunov 2007 eq. 19
    const Real A1 = 2.7822_rt;
    const Real A2 = 98.34_rt;
    const Real A3 = sqrt(3.0_rt) - A1 / sqrt(A2);
    const Real B1 = -1.7476_rt;
    const Real B2 = 66.07_rt;
    const Real B3 = 1.12_rt;
    const Real B4 = 65_rt;
    number_t gamtilde2 = gamtilde * gamtilde;

    // Chugunov 2007 eq. 19
    number_t term1 = 1.0_rt / sqrt(A2 + gamtilde);
    number_t term2 = 1.0_rt / (1.0_rt + gamtilde);
    number_t term3 = gamtilde2 / (B2 + gamtilde);
    number_t term4 = gamtilde2 / (B4 + gamtilde2);

    number_t inner = A1 * term1 + A3 * term2;
    number_t h = pow(gamtilde, 1.5_rt) * inner + B1 * term3 + B3 * term4;

    // Assume h >= 0
    return max(h, 0.0_rt);
}


template <typename number_t>
inline
number_t actual_log_screen(const plasma_state_t<number_t>& state,
                         const scrn::screen_factors_t& scn_fac)
{
    return chugunov2007(state, scn_fac);
}

inline
void actual_log_screen(const plasma_state_t<Real>& state,
                       const scrn::screen_factors_t& scn_fac,
                       Real& log_scor)
{
    log_scor = chugunov2007(state, scn_fac);
}

// Expose both composition partials from one screening evaluation.
inline
void actual_log_screen(const plasma_state_t<screening_dual_t>& state,
                       const scrn::screen_factors_t& scn_fac,
                       Real& log_scor, Real& dlog_scor_dYe, Real& dlog_scor_dS2)
{
    const auto result = actual_log_screen(state, scn_fac);
    log_scor = result.value;
    dlog_scor_dYe = result.deriv.d1;
    dlog_scor_dS2 = result.deriv.d2;
}

template <typename number_t>
inline
number_t actual_screen(const plasma_state_t<number_t>& state,
                       const scrn::screen_factors_t& scn_fac)
{
    using std::exp;
    using std::min;

    number_t log_scor = actual_log_screen(state, scn_fac);
    if (log_scor <= 0.0_rt) {
        return 1.0_rt;
    }
    // Machine limit the output, including zero derivatives when capped.
    constexpr Real log_scor_max = 300.e0_rt;
    return exp(min(log_scor, log_scor_max));
}

inline
void actual_screen(const plasma_state_t<Real>& state,
                   const scrn::screen_factors_t& scn_fac,
                   Real& scor)
{
    scor = actual_screen(state, scn_fac);
}

#endif
