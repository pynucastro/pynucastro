#include <actual_network.H>


#ifdef NSE_NET
namespace NSE_INDEX
{
    AMREX_GPU_MANAGED amrex::Array2D<std::int8_t, 1, Rates::NumRates, 1, 7, amrex::Order::C> rate_indices {
        -1, -1, -1, -1, -1, -1, -1,  // He4_Mg24_to_Si28_removed
        -1, -1, -1, -1, -1, -1, -1,  // He4_Mg24_to_p_Al27_removed
        -1, -1, -1, -1, -1, -1, -1,  // p_Al27_to_Si28_removed
        -1, -1, -1, -1, -1, -1, -1,  // Si28_to_He4_Mg24_removed
        -1, -1, -1, -1, -1, -1, -1,  // Si28_to_p_Al27_removed
        -1, -1, -1, -1, -1, -1, -1,  // p_Al27_to_He4_Mg24_removed
        -1, -1, -1, -1, -1, -1, -1,  // He4_Si28_to_S32_removed
        -1, -1, -1, -1, -1, -1, -1,  // He4_Si28_to_p_P31_removed
        -1, -1, -1, -1, -1, -1, -1,  // p_P31_to_S32_removed
        -1, -1, -1, -1, -1, -1, -1,  // S32_to_He4_Si28_removed
        -1, -1, -1, -1, -1, -1, -1,  // S32_to_p_P31_removed
        -1, -1, -1, -1, -1, -1, -1,  // p_P31_to_He4_Si28_removed
        -1, 0, 1, -1, -1, 2, 14,  // Mg24_He4_to_Si28_approx
        -1, -1, 2, -1, 0, 1, -1,  // Si28_to_Mg24_He4_approx
        -1, 0, 2, -1, -1, 3, 16,  // Si28_He4_to_S32_approx
        -1, -1, 3, -1, 0, 2, -1  // S32_to_Si28_He4_approx
    };
}
#endif

void actual_network_init()
{

}
