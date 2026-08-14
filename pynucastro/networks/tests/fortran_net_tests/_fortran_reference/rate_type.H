#ifndef NET_RATE_TYPE_H
#define NET_RATE_TYPE_H

#include <actual_network.H>

struct rate_t {
    Array1D<Real, 1, Rates::NumRates>  screened_rates;
#ifdef SCREENING
    Array1D<Real, 1, Rates::NumScreenPairs>  log_screen;
#endif
    Real enuc_weak;
};

struct rate_derivs_t {
    Array1D<Real, 1, Rates::NumRates>  screened_rates;
    Array1D<Real, 1, Rates::NumRates>  dscreened_rates_dT;
#ifdef SCREENING
    Array1D<Real, 1, Rates::NumScreenPairs>  log_screen;
    Array1D<Real, 1, Rates::NumScreenPairs>  dlog_screen_dT;
#endif
    Real enuc_weak;
};


#endif