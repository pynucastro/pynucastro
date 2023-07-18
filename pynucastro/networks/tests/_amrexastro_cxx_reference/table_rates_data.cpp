#include <AMReX_Array.H>
#include <string>
#include <table_rates.H>
#include <AMReX_Print.H>

using namespace amrex;

namespace rate_tables
{

    AMREX_GPU_MANAGED table_t j_na23_ne23_meta;
    AMREX_GPU_MANAGED Array3D<Real, 1, 39, 1, 152, 1, 6> j_na23_ne23_data;
    AMREX_GPU_MANAGED Array1D<Real, 1, 152> j_na23_ne23_rhoy;
    AMREX_GPU_MANAGED Array1D<Real, 1, 39> j_na23_ne23_temp;

    AMREX_GPU_MANAGED table_t j_ne23_na23_meta;
    AMREX_GPU_MANAGED Array3D<Real, 1, 39, 1, 152, 1, 6> j_ne23_na23_data;
    AMREX_GPU_MANAGED Array1D<Real, 1, 152> j_ne23_na23_rhoy;
    AMREX_GPU_MANAGED Array1D<Real, 1, 39> j_ne23_na23_temp;


}


void init_tabular()
{

    amrex::Print() << "reading in network electron-capture / beta-decay tables..." << std::endl;

    using namespace rate_tables;

    j_na23_ne23_meta.ntemp = 39;
    j_na23_ne23_meta.nrhoy = 152;
    j_na23_ne23_meta.nvars = 6;
    j_na23_ne23_meta.nheader = 7;

    init_tab_info(j_na23_ne23_meta, "23na-23ne_electroncapture.dat", j_na23_ne23_rhoy, j_na23_ne23_temp, j_na23_ne23_data);


    j_ne23_na23_meta.ntemp = 39;
    j_ne23_na23_meta.nrhoy = 152;
    j_ne23_na23_meta.nvars = 6;
    j_ne23_na23_meta.nheader = 5;

    init_tab_info(j_ne23_na23_meta, "23ne-23na_betadecay.dat", j_ne23_na23_rhoy, j_ne23_na23_temp, j_ne23_na23_data);



}
