#ifndef screen_data_H
#define screen_data_H

#include <amrex_bridge.H>
#include <network_properties.H>

namespace scrn {
    class screen_factors_t {

    public:

        Real z1 = -1;
        Real z2 = -1;
        Real a1 = -1;
        Real a2 = -1;

        // ztilde  = effective ion radius factor for a MCP
        // ztilde3 = ztilde**3

        Real ztilde = 0.0;
        Real ztilde3 = 0.0;

    };

    inline
    screen_factors_t calculate_screen_factor (Real z1, Real a1,
                                              Real z2, Real a2)
    {
        screen_factors_t scn_fac;

        scn_fac.z1 = z1;
        scn_fac.a1 = a1;
        scn_fac.z2 = z2;
        scn_fac.a2 = a2;

        scn_fac.ztilde = 0.5_rt * (std::pow(z1, 1.0_rt / 3.0_rt) +
                                   std::pow(z2, 1.0_rt / 3.0_rt));
        scn_fac.ztilde3 = std::pow(scn_fac.ztilde, 3.0);

        return scn_fac;
    }
}

#endif
