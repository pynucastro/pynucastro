"""A C++ reaction network for integration into the AMReX Astro
Microphysics set of reaction networks used by astrophysical hydrodynamics
codes"""


import glob
import os

from pynucastro.networks.base_cxx_network import BaseCxxNetwork
from pynucastro.nucdata import Nucleus
from pynucastro.rates import ReacLibRate


class AmrexAstroCxxNetwork(BaseCxxNetwork):
    def __init__(self, *args, **kwargs):

        # this network can have a special kwarg called disable_rate_params
        try:
            disable_rate_params = kwargs.pop("disable_rate_params")
        except KeyError:
            disable_rate_params = []

        # Initialize BaseCxxNetwork parent class
        super().__init__(*args, **kwargs)

        self.ftags['<rate_param_tests>'] = self._rate_param_tests
        self.ftags['<rate_indices>'] = self._fill_rate_indices
        self.ftags['<npa_index>'] = self._fill_npa_index

        self.disable_rate_params = disable_rate_params
        self.function_specifier = "AMREX_GPU_HOST_DEVICE AMREX_INLINE"
        self.dtype = "Real"

    def _get_template_files(self):

        template_pattern = os.path.join(self.pynucastro_dir,
                                        'templates',
                                        'amrexastro-cxx-microphysics',
                                        '*.template')

        return glob.glob(template_pattern)

    def _rate_param_tests(self, n_indent, of):

        for _, r in enumerate(self.rates):
            if r in self.disable_rate_params:
                of.write(f"{self.indent*n_indent}if (disable_{r.cname()}) {{\n")
                of.write(f"{self.indent*n_indent}    rate_eval.screened_rates(k_{r.cname()}) = 0.0;\n")
                of.write(f"{self.indent*n_indent}    if constexpr (std::is_same<T, rate_derivs_t>::value) {{\n")
                of.write(f"{self.indent*n_indent}        rate_eval.dscreened_rates_dT(k_{r.cname()}) = 0.0;\n")
                of.write(f"{self.indent*n_indent}    }}\n")
                # check for the reverse too -- we disable it with the same parameter
                rr = self.find_reverse(r)
                if rr is not None:
                    of.write(f"{self.indent*n_indent}    rate_eval.screened_rates(k_{rr.cname()}) = 0.0;\n")
                    of.write(f"{self.indent*n_indent}    if constexpr (std::is_same<T, rate_derivs_t>::value) {{\n")
                    of.write(f"{self.indent*n_indent}        rate_eval.dscreened_rates_dT(k_{rr.cname()}) = 0.0;\n")
                    of.write(f"{self.indent*n_indent}    }}\n")
                of.write(f"{self.indent*n_indent}}}\n\n")

    def _write_network(self, odir=None):
        """
        This writes the RHS, jacobian and ancillary files for the system of ODEs that
        this network describes, using the template files.
        """

        super()._write_network(odir=odir)

        if odir is None:
            odir = os.getcwd()
        # create a .net file with the nuclei properties
        with open(os.path.join(odir, "pynucastro.net"), "w") as of:
            for nuc in self.unique_nuclei:
                of.write(f"{nuc.spec_name:25} {nuc.short_spec_name:6} {nuc.A:6.1f} {nuc.Z:6.1f}\n")

            for nuc in self.approx_nuclei:
                of.write(f"__extra_{nuc.spec_name:17} {nuc.short_spec_name:6} {nuc.A:6.1f} {nuc.Z:6.1f}\n")

        # write the _parameters file
        with open(os.path.join(odir, "_parameters"), "w") as of:
            of.write("@namespace: network\n\n")
            if self.disable_rate_params:
                for r in self.disable_rate_params:
                    of.write(f"disable_{r.cname()}    int     0\n")

    def _fill_npa_index(self, n_indent, of):
        #Get the index of h1, neutron, and helium-4 if they're present in the network.

        LIG = list(map(Nucleus, ["p", "n", "he4"]))

        for nuc in LIG:
            if nuc in self.unique_nuclei:
                of.write(f"{self.indent*n_indent}constexpr int {nuc.short_spec_name}_index = {self.unique_nuclei.index(nuc)};\n")
            else:
                of.write(f"{self.indent*n_indent}constexpr int {nuc.short_spec_name}_index = -1;\n")

    def _fill_rate_indices(self, n_indent, of):
        """
        Fills the index needed for the NSE_NET algorithm.

        Fill rate_indices: 2D array with 1-based index of shape of size (NumRates, 7).
           - Each row represents a rate in self.all_rates.
           - The first 3 elements of the row represents the index of reactants in self.unique_nuclei
           - The next 3 elements of the row represents the index of the products in self.unique_nuclei.
           - The 7th element of the row represents the index of the corresponding reverse rate
             (set to -1 if no corresponding reverse rate). This is a 1-based instead of 0-based index.
           - Set all elements of the current row to -1 if the rate has removed suffix
             indicating its not directly in the network.
        """

        # Fill in the rate indices
        of.write(f"{self.indent*n_indent}AMREX_GPU_MANAGED amrex::Array2D<int, 1, Rates::NumRates, 1, 7, Order::C> rate_indices {{\n")

        for n, rate in enumerate(self.all_rates):
            tmp = ','
            if n == len(self.all_rates) - 1:
                tmp = ''

            # meaning it is removed.
            if isinstance(rate, ReacLibRate) and rate.removed is not None:
                of.write(f"{self.indent*n_indent}    -1, -1, -1, -1, -1, -1, -1{tmp}\n")
                continue

            # Find the reactants and products indices
            reactant_ind = [-1 for n in range(3 - len(rate.reactants))]
            product_ind = [-1 for n in range(3 - len(rate.products))]

            for nuc in rate.reactants:
                reactant_ind.append(self.unique_nuclei.index(nuc))

            for nuc in rate.products:
                product_ind.append(self.unique_nuclei.index(nuc))

            reactant_ind.sort()
            product_ind.sort()

            # Find the reverse rate index
            rr_ind = -1
            rr = self.find_reverse(rate)

            # Note that rate index is 1-based
            if rr is not None:
                rr_ind = self.all_rates.index(rr) + 1

            of.write(f"{self.indent*n_indent}    {reactant_ind[0]}, {reactant_ind[1]}, {reactant_ind[2]}, {product_ind[0]}, {product_ind[1]}, {product_ind[2]}, {rr_ind}{tmp}\n")

        of.write(f"{self.indent*n_indent}}};\n")
