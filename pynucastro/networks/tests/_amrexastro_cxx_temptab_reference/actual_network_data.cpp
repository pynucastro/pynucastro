#include <actual_network.H>


#ifdef NSE_NET
namespace NSE_INDEX
{
    AMREX_GPU_MANAGED amrex::Array2D<std::int8_t, 1, Rates::NumRates, 1, 7, amrex::Order::C> rate_indices {
        -1, 0, 1, -1, -1, 2, 2,  // O16_p_to_F17_iliadis
        -1, -1, 2, -1, 0, 1, -1  // F17_to_O16_p_derived
    };
}
#endif

void actual_network_init()
{

}
