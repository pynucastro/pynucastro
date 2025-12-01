#include <actual_network.H>


#ifdef NSE_NET
namespace NSE_INDEX
{
    AMREX_GPU_MANAGED amrex::Array2D<std::uint8_t, 1, Rates::NumRates, 1, 7, amrex::Order::C> rate_indices {
        -1, 1, 2, -1, -1, 4, 4,  // He4_Fe52__Ni56
        -1, 0, 3, -1, -1, 4, 5,  // p_Co55__Ni56
        -1, 1, 2, -1, 0, 3, 6,  // He4_Fe52__p_Co55
        -1, -1, 4, -1, 1, 2, -1,  // Ni56__He4_Fe52__derived
        -1, -1, 4, -1, 0, 3, -1,  // Ni56__p_Co55__derived
        -1, 0, 3, -1, 1, 2, -1  // p_Co55__He4_Fe52__derived
    };
}
#endif

void actual_network_init()
{

}
