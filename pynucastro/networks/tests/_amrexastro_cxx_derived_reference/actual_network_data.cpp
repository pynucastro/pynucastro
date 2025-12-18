#include <actual_network.H>


#ifdef NSE_NET
namespace NSE_INDEX
{
    AMREX_GPU_MANAGED amrex::Array2D<std::int8_t, 1, Rates::NumRates, 1, 7, amrex::Order::C> rate_indices {
        -1, 1, 2, -1, -1, 4, 4,  // He4_Fe52_to_Ni56_reaclib
        -1, 0, 3, -1, -1, 4, 5,  // p_Co55_to_Ni56_reaclib
        -1, 1, 2, -1, 0, 3, 6,  // He4_Fe52_to_p_Co55_reaclib
        -1, -1, 4, -1, 1, 2, -1,  // Ni56_to_He4_Fe52_derived
        -1, -1, 4, -1, 0, 3, -1,  // Ni56_to_p_Co55_derived
        -1, 0, 3, -1, 1, 2, -1  // p_Co55_to_He4_Fe52_derived
    };
}
#endif

void actual_network_init()
{

}
