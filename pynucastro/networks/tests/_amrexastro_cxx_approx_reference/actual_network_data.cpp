#include <actual_network.H>


#ifdef NSE_NET
namespace NSE_INDEX
{
    AMREX_GPU_MANAGED amrex::Array2D<std::uint8_t, 1, Rates::NumRates, 1, 7, amrex::Order::C> rate_indices {
        -1, -1, -1, -1, -1, -1, -1,  // He4_Mg24__Si28__removed
        -1, -1, -1, -1, -1, -1, -1,  // He4_Mg24__p_Al27__removed
        -1, -1, -1, -1, -1, -1, -1,  // p_Al27__Si28__removed
        -1, -1, -1, -1, -1, -1, -1,  // Si28__He4_Mg24__removed
        -1, -1, -1, -1, -1, -1, -1,  // Si28__p_Al27__removed
        -1, -1, -1, -1, -1, -1, -1,  // p_Al27__He4_Mg24__removed
        -1, -1, -1, -1, -1, -1, -1,  // He4_Si28__S32__removed
        -1, -1, -1, -1, -1, -1, -1,  // He4_Si28__p_P31__removed
        -1, -1, -1, -1, -1, -1, -1,  // p_P31__S32__removed
        -1, -1, -1, -1, -1, -1, -1,  // S32__He4_Si28__removed
        -1, -1, -1, -1, -1, -1, -1,  // S32__p_P31__removed
        -1, -1, -1, -1, -1, -1, -1,  // p_P31__He4_Si28__removed
        -1, 0, 1, -1, -1, 2, 14,  // Mg24_He4__Si28__approx
        -1, -1, 2, -1, 0, 1, -1,  // Si28__Mg24_He4__approx
        -1, 0, 2, -1, -1, 3, 16,  // Si28_He4__S32__approx
        -1, -1, 3, -1, 0, 2, -1  // S32__Si28_He4__approx
    };
}
#endif

void actual_network_init()
{

}
