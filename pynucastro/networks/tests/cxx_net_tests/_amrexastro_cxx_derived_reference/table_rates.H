#ifndef TABLE_RATES_H
#define TABLE_RATES_H

#include <iostream>
#include <sstream>
#include <fstream>
#include <string>

#include <AMReX_Array.H>
#include <AMReX_Math.H>
#include <interp_tools.H>

using namespace amrex::literals;


// Table is expected to be in terms of dens*ye and temp (logarithmic, cgs units)
// Table energy units are expected in terms of ergs

// all tables are expected to have columns:
// Log(rhoY)     Log(T)   mu    dQ    Vs    Log(e-cap-rate)   Log(nu-energy-loss)  Log(gamma-energy)
// Log(g/cm^3)   Log(K)   erg   erg   erg   Log(1/s)          Log(erg/s)           Log(erg/s)
//

constexpr int num_tables = 0;

enum TableVars : std::uint8_t
{
    jtab_rate    = 1,
    jtab_nuloss  = 2,
    jtab_gamma   = 3,
    num_vars = jtab_gamma
};


struct table_t
{
    int ntemp;
    int nrhoy;
    int nvars;
    int nheader;
};

// we add two additional indices:
//    k_index_dlogr_dlogT : d Log(rate) / d Log(temperature)
//    k_index_dlogenuc_dlogT : d Log(nuloss) / d Log(temperature)
//
// These are not indices nto the table but into the 'entries'
// array. Is important to mention that although we compute dlogr/dlogT
// as the computed quantity in 'entries', we pursue ultimately dr/dt
// as the final desired quantity to be computed for this index.

constexpr int k_index_dlogr_dlogT = num_vars + 1;
constexpr int k_index_dlognu_dlogT = num_vars + 2;
constexpr int add_vars = 2;  // 1 Additional Var in entries


namespace rate_tables
{
}


AMREX_INLINE AMREX_GPU_DEVICE
amrex::Real
evaluate_linear_1d(const amrex::Real fhi, const amrex::Real flo, const amrex::Real xhi, const amrex::Real xlo, const amrex::Real x)
{
    // This function is a 1-D linear interpolator, that keeps x constant to xlo or xhi, based
    // on the side, if x is outside [xlo, xhi] respectively.

    amrex::Real xx = amrex::Clamp(x, xlo, xhi);
    amrex::Real f = flo + (fhi - flo) * (xx - xlo) / (xhi - xlo);

    return f;
}

AMREX_INLINE AMREX_GPU_DEVICE
amrex::Real
evaluate_linear_2d(const amrex::Real fip1jp1, const amrex::Real fip1j, const amrex::Real fijp1, const amrex::Real fij,
          const amrex::Real xhi, const amrex::Real xlo, const amrex::Real yhi, const amrex::Real ylo,
          const amrex::Real x, const amrex::Real y)
{
    // This is the 2-D linear interpolator, as an extension of evaluate_linear_1d.

    const amrex::Real u = (amrex::Clamp(x, xlo, xhi) - xlo) / (xhi - xlo);
    const amrex::Real v = (amrex::Clamp(y, ylo, yhi) - ylo) / (yhi - ylo);

    const amrex::Real E = fij;
    const amrex::Real C = fijp1 - fij;
    const amrex::Real B = fip1j - fij;
    const amrex::Real A = fip1jp1 - fip1j - fijp1 + fij;

    return A * u * v + B * u + C * v + E;
}


template<typename R, typename T, typename D>
AMREX_INLINE AMREX_GPU_DEVICE
amrex::Real
evaluate_vars(const int irhoy_lo, const int jtemp_lo,
              const R& log_rhoy_table, const T& log_temp_table, const D& data,
              const amrex::Real log_rhoy, const amrex::Real log_temp, const int component)
{
    // This function evaluates the 2-D interpolator, for several pairs of rho_ye and temperature.

    int irhoy_hi = irhoy_lo + 1;
    int jtemp_hi = jtemp_lo + 1;

    amrex::Real rhoy_lo = log_rhoy_table(irhoy_lo);
    amrex::Real rhoy_hi = log_rhoy_table(irhoy_hi);

    amrex::Real t_lo = log_temp_table(jtemp_lo);
    amrex::Real t_hi = log_temp_table(jtemp_hi);

    amrex::Real fij     = data(jtemp_lo, irhoy_lo, component);
    amrex::Real fip1j   = data(jtemp_lo, irhoy_hi, component);
    amrex::Real fijp1   = data(jtemp_hi, irhoy_lo, component);
    amrex::Real fip1jp1 = data(jtemp_hi, irhoy_hi, component);

    amrex::Real r = evaluate_linear_2d(fip1jp1, fip1j, fijp1, fij,
                                       rhoy_hi, rhoy_lo, t_hi, t_lo, log_rhoy, log_temp);

    return r;
}


template<typename R, typename T, typename D>
AMREX_INLINE AMREX_GPU_DEVICE
amrex::Real
evaluate_dr_dtemp(const int irhoy_lo, const int jtemp_lo,
                  const table_t& table_meta, const R& log_rhoy_table, const T& log_temp_table, const D& data,
                  const amrex::Real log_rhoy, const amrex::Real log_temp, const int component)
{

    // Differentiate the same bilinear log-rate interpolant as evaluate_vars.
    // At temperature knots, use the cell selected by find_index_extrap:
    // the cell to the right, except at the upper endpoint (the last cell).

    const int irhoy_hi = irhoy_lo + 1;
    const int jtemp_hi = jtemp_lo + 1;

    const amrex::Real log_temp_lo = log_temp_table(jtemp_lo);
    const amrex::Real log_temp_hi = log_temp_table(jtemp_hi);

    // The interpolated rate is constant outside the temperature table.
    if (log_temp < log_temp_lo || log_temp > log_temp_hi) {
        return 0.0_rt;
    }

    const amrex::Real log_rhoy_lo = log_rhoy_table(irhoy_lo);
    const amrex::Real log_rhoy_hi = log_rhoy_table(irhoy_hi);

    const amrex::Real fij     = data(jtemp_lo, irhoy_lo, component);
    const amrex::Real fip1j   = data(jtemp_lo, irhoy_hi, component);
    const amrex::Real fijp1   = data(jtemp_hi, irhoy_lo, component);
    const amrex::Real fip1jp1 = data(jtemp_hi, irhoy_hi, component);

    // compute the gradients in the T direction: dr/dT |_i and dr/dT |_{i+1}

    const amrex::Real dlogT = log_temp_hi - log_temp_lo;
    const amrex::Real dlogr_dlogT_i = (fijp1 - fij) / dlogT;
    const amrex::Real dlogr_dlogT_ip1 = (fip1jp1 - fip1j) / dlogT;

    return evaluate_linear_1d(dlogr_dlogT_ip1, dlogr_dlogT_i,
                              log_rhoy_hi, log_rhoy_lo, log_rhoy);
}


template <int do_T_derivatives, typename R, typename T, typename D>
AMREX_INLINE AMREX_GPU_DEVICE
void
get_entries(const int irhoy_lo, const int jtemp_lo,
            const table_t& table_meta, const R& log_rhoy_table, const T& log_temp_table, const D& data,
            const amrex::Real log_rhoy, const amrex::Real log_temp, amrex::Array1D<amrex::Real, 1, num_vars+add_vars>& entries)
{

    // only interpolate the quantities we currently use

    entries(jtab_rate) = evaluate_vars(irhoy_lo, jtemp_lo,
                                       log_rhoy_table, log_temp_table, data,
                                       log_rhoy, log_temp, jtab_rate);

    entries(jtab_nuloss) = evaluate_vars(irhoy_lo, jtemp_lo,
                                         log_rhoy_table, log_temp_table, data,
                                         log_rhoy, log_temp, jtab_nuloss);

    entries(jtab_gamma) = evaluate_vars(irhoy_lo, jtemp_lo,
                                        log_rhoy_table, log_temp_table, data,
                                        log_rhoy, log_temp, jtab_gamma);

    if constexpr (do_T_derivatives) {
        entries(k_index_dlogr_dlogT) =
            evaluate_dr_dtemp(irhoy_lo, jtemp_lo,
                              table_meta, log_rhoy_table, log_temp_table, data,
                              log_rhoy, log_temp, jtab_rate);

        entries(k_index_dlognu_dlogT) =
            evaluate_dr_dtemp(irhoy_lo, jtemp_lo,
                              table_meta, log_rhoy_table, log_temp_table, data,
                              log_rhoy, log_temp, jtab_nuloss);
    }
}

template <int do_T_derivatives, typename R, typename T, typename D>
AMREX_INLINE AMREX_GPU_DEVICE
void
tabular_evaluate(const table_t& table_meta,
                 const R& log_rhoy_table, const T& log_temp_table, const D& data,
                 const amrex::Real log_rhoy, const amrex::Real log_temp, const amrex::Real temp,
                 amrex::Real& rate, amrex::Real& drate_dT,
                 amrex::Real& enu, amrex::Real& denu_dT, amrex::Real& edot_gamma)
{
    amrex::Array1D<amrex::Real, 1, num_vars+add_vars> entries;

    // Get the table entries at this rhoy, temp

    // get the index into the table -- all components are
    // stored with the same indices

    const int irhoy_lo = interp_net::find_index_extrap(log_rhoy, log_rhoy_table);
    const int jtemp_lo = interp_net::find_index_extrap(log_temp, log_temp_table);

    get_entries<do_T_derivatives>(irhoy_lo, jtemp_lo,
                                  table_meta, log_rhoy_table, log_temp_table, data,
                                  log_rhoy, log_temp, entries);

    // Fill outputs: rate, d(rate)/d(temperature), and
    // (negative) neutrino loss contribution to energy generation

    rate = amrex::Math::exp10(entries(jtab_rate));
    if constexpr (do_T_derivatives) {
        drate_dT = rate * entries(k_index_dlogr_dlogT) / temp;
    } else {
        drate_dT = 0.0;
    }
    enu = -amrex::Math::exp10(entries(jtab_nuloss));
    if constexpr (do_T_derivatives) {
        denu_dT = enu * entries(k_index_dlognu_dlogT) / temp;
    } else {
        denu_dT = 0.0;
    }
    edot_gamma = amrex::Math::exp10(entries(jtab_gamma));
}

#endif
