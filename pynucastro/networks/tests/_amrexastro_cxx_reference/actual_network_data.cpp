#include <actual_network.H>


#ifdef NSE_NET
namespace NSE_INDEX
{
    AMREX_GPU_MANAGED amrex::Array2D<int, 1, Rates::NumRates, 1, 7, amrex::Order::C> rate_indices {
        -1, 3, 3, -1, 2, 5, -1,
        -1, 3, 3, -1, 0, 8, -1,
        -1, 3, 3, -1, 1, 7, -1,
        -1, 2, 3, -1, -1, 4, -1,
        -1, -1, 0, -1, -1, 1, -1,
        -1, -1, 7, -1, -1, 6, -1,
        -1, -1, 6, -1, -1, 7, 6
    };
}
#endif

void actual_network_init()
{

}
