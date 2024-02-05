#include <actual_network.H>


#ifdef NSE_NET
namespace NSE_INDEX
{
    AMREX_GPU_MANAGED amrex::Array2D<int, 1, Rates::NumRates, 1, 7, Order::C> rate_indices {
        -1, 1, 2, -1, -1, 4, 4,
        -1, 0, 3, -1, -1, 4, 5,
        -1, 1, 2, -1, 0, 3, 6,
        -1, -1, 4, -1, 1, 2, -1,
        -1, -1, 4, -1, 0, 3, -1,
        -1, 0, 3, -1, 1, 2, -1
    };
}
#endif

void actual_network_init()
{

}
