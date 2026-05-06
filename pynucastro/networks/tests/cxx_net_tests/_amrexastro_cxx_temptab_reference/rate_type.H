#ifndef NET_RATE_TYPE_H
#define NET_RATE_TYPE_H

#include <actual_network.H>

struct rate_t {
    amrex::Array1D<amrex::Real, 1, Rates::NumRates>  screened_rates;
#ifdef SCREENING
    amrex::Array1D<amrex::Real, 1, Rates::NumScreenPairs>  log_screen;
#endif
    amrex::Real enuc_weak;
};

struct rate_derivs_t {
    amrex::Array1D<amrex::Real, 1, Rates::NumRates>  screened_rates;
    amrex::Array1D<amrex::Real, 1, Rates::NumRates>  dscreened_rates_dT;
#ifdef SCREENING
    amrex::Array1D<amrex::Real, 1, Rates::NumScreenPairs>  log_screen;
    amrex::Array1D<amrex::Real, 1, Rates::NumScreenPairs>  dlog_screen_dT;
#endif
    amrex::Real enuc_weak;
};


#endif