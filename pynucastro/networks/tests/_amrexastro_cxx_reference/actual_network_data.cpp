#include <actual_network.H>


#ifdef NSE_NET
namespace NSE_INDEX
{
    AMREX_GPU_MANAGED amrex::Array2D<std::int8_t, 1, Rates::NumRates, 1, 7, amrex::Order::C> rate_indices {
        -1, 3, 3, -1, 2, 5, -1,  // C12_C12_to_He4_Ne20
        -1, 3, 3, -1, 0, 8, -1,  // C12_C12_to_n_Mg23
        -1, 3, 3, -1, 1, 7, -1,  // C12_C12_to_p_Na23
        -1, 2, 3, -1, -1, 4, -1,  // He4_C12_to_O16
        -1, -1, 0, -1, -1, 1, -1,  // n_to_p_weak_wc12
        -1, -1, 7, -1, -1, 6, -1,  // Na23_to_Ne23
        -1, -1, 6, -1, -1, 7, 6  // Ne23_to_Na23
    };
}
#endif

void actual_network_init()
{

}
