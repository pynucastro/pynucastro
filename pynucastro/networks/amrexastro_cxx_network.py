"""A C++ reaction network for integration into the AMReX Astro
Microphysics set of reaction networks used by astrophysical
hydrodynamics codes

"""

import itertools
import re
from pathlib import Path

from pynucastro.constants import constants
from pynucastro.networks.base_cxx_network import (BaseCxxNetwork,
                                                  _signed_rate_dtype)
from pynucastro.nucdata import Nucleus


class AmrexAstroCxxNetwork(BaseCxxNetwork):
    """A C++ network for simulation codes based on AMReX."""

    def __init__(self, *args, **kwargs):

        # this network can have a special kwarg called disable_rate_params
        try:
            disable_rate_params = kwargs.pop("disable_rate_params")
        except KeyError:
            disable_rate_params = []

        # Initialize BaseCxxNetwork parent class
        super().__init__(*args, **kwargs)

        self.ftags['<rate_param_tests>'] = self._rate_param_tests
        self.ftags['<npa_index>'] = self._fill_npa_index
        self.ftags['<nse_rate_pair_data>'] = self._write_nse_rate_pair_data

        self.disable_rate_params = disable_rate_params
        self.function_specifier = "AMREX_GPU_HOST_DEVICE AMREX_INLINE"
        self.dtype = "amrex::Real"
        self.array_namespace = "amrex::"

    def _get_template_files(self):

        path = self.pynucastro_dir/"templates/amrexastro-cxx-microphysics"

        return path.glob("*.template")

    def _rate_param_tests(self, n_indent, of):

        for _, r in enumerate(self.rates):
            if r in self.disable_rate_params:
                of.write(f"{self.indent*n_indent}if (disable_{r.fname}) {{\n")
                of.write(f"{self.indent*n_indent}    rate_eval.screened_rates(k_{r.fname}) = 0.0;\n")
                of.write(f"{self.indent*n_indent}    if constexpr (std::is_same_v<T, rate_derivs_t>) {{\n")
                of.write(f"{self.indent*n_indent}        rate_eval.dscreened_rates_dT(k_{r.fname}) = 0.0;\n")
                of.write(f"{self.indent*n_indent}    }}\n")
                # check for the reverse too -- we disable it with the same parameter
                rr = self.find_reverse(r)
                if rr is not None:
                    of.write(f"{self.indent*n_indent}    rate_eval.screened_rates(k_{rr.fname}) = 0.0;\n")
                    of.write(f"{self.indent*n_indent}    if constexpr (std::is_same_v<T, rate_derivs_t>) {{\n")
                    of.write(f"{self.indent*n_indent}        rate_eval.dscreened_rates_dT(k_{rr.fname}) = 0.0;\n")
                    of.write(f"{self.indent*n_indent}    }}\n")
                of.write(f"{self.indent*n_indent}}}\n\n")

    def _ebind(self, n_indent, of):
        for n, nuc in enumerate(self.unique_nuclei):
            if n == 0:
                of.write(f"{self.indent*n_indent}if constexpr (spec == {nuc.cindex()}) {{\n")
            else:
                of.write(f"{self.indent*n_indent}else if constexpr (spec == {nuc.cindex()}) {{\n")
            of.write(f"{self.indent*(n_indent+1)}return {nuc.nucbind * nuc.A}_rt;\n")
            of.write(f"{self.indent*(n_indent)}}}\n")

    def _mion(self, n_indent, of):
        for n, nuc in enumerate(self.unique_nuclei):
            if n == 0:
                of.write(f"{self.indent*n_indent}if constexpr (spec == {nuc.cindex()}) {{\n")
            else:
                of.write(f"{self.indent*n_indent}else if constexpr (spec == {nuc.cindex()}) {{\n")
            of.write(f"{self.indent*(n_indent+1)}return {nuc.A_nuc * constants.m_u_C18}_rt;\n")
            of.write(f"{self.indent*(n_indent)}}}\n")

    def _fill_spin_state_cases(self, n_indent, of):
        def key_func(nuc):
            if nuc.spin_states is None:
                return -1
            return nuc.spin_states

        # group identical cases together to satisfy clang-tidy
        nuclei = sorted(self.unique_nuclei + self.approx_nuclei, key=key_func)

        FIRST_ENCOUNTER = True
        for spin_state, group in itertools.groupby(nuclei, key=key_func):
            if spin_state == -1:
                continue

            if FIRST_ENCOUNTER:
                of.write(f"{self.indent*n_indent}if constexpr (\n")
                parenthesis_indent = f"{self.indent*(n_indent+1)}          "
                FIRST_ENCOUNTER = False
            else:
                of.write(f"{self.indent*n_indent}else if constexpr (\n")
                parenthesis_indent = f"{self.indent*(n_indent+1)}               "

            # Divide group of spec into subgroups of 3 for better formatting
            group = list(group)
            subgroups = [group[n:n+3] for n in range(0, len(group), 3)]
            for i, subgroup in enumerate(subgroups):
                spec_string = " || ".join([f"spec == {n.cindex()}" for n in subgroup])

                # If it is not the last subgroup, add || in the end
                if i != len(subgroups) - 1:
                    spec_string += " ||"
                of.write(f"{self.indent*(n_indent+1)}{spec_string}\n")

            of.write(f"{parenthesis_indent})\n")
            of.write(f"{self.indent*n_indent}{{\n")
            of.write(f"{self.indent*(n_indent+1)}return {spin_state}.0_rt;\n")
            of.write(f"{self.indent*n_indent}}}\n")

    def _cxxify(self, s):
        # Replace std::pow(x, n) with amrex::Math::powi<n>(x) for amrexastro_cxx_network

        cxx_code = super()._cxxify(s)
        std_pow_pattern = r"std::pow\(([^,]+),\s*(\d+)\)"
        amrex_powi_replacement = r"amrex::Math::powi<\2>(\1)"
        return re.sub(std_pow_pattern, amrex_powi_replacement, cxx_code)

    def _write_network(self, odir=None):
        """Output the RHS, jacobian and ancillary files for the system
        of ODEs that this network describes, using the template files.

        """

        super()._write_network(odir=odir)

        if odir is None:
            odir = Path.cwd()
        # create a .net file with the nuclei properties
        with open(Path(odir, "pynucastro.net"), "w") as of:
            for nuc in self.unique_nuclei:
                short_spec_name = nuc.short_spec_name
                if nuc.short_spec_name != "n":
                    short_spec_name = nuc.short_spec_name.capitalize()
                of.write(f"{nuc.spec_name:25} {short_spec_name:6} {nuc.A:6.1f} {nuc.Z:6.1f}\n")

            for nuc in self.approx_nuclei:
                short_spec_name = nuc.short_spec_name
                if nuc.short_spec_name != "n":
                    short_spec_name = nuc.short_spec_name.capitalize()
                of.write(f"__extra_{nuc.spec_name:17} {short_spec_name:6} {nuc.A:6.1f} {nuc.Z:6.1f}\n")

        # write the _parameters file
        with open(Path(odir, "_parameters"), "w") as of:
            of.write("@namespace: network\n\n")
            if self.disable_rate_params:
                for r in self.disable_rate_params:
                    of.write(f"disable_{r.fname}    int     0\n")

    def _get_nse_rate_pairs(self):
        """Return a list of RatePairs eligible for NSE_NET grouping.

        Filtering criteria:
        1. Must have a reverse rate (RatePairs only)
        2. Excludes weak or removed rates
        3. At most 2 reactants and 2 products (except for triple-alpha)
        4. Only 1 or 2 nuclei not in {p, n, He4} across reactants + products

        """

        nse_rate_pairs = []

        # Light isotope group
        LIG = Nucleus.cast_list(["p", "n", "he4"])

        # Triple alpha reactants and products
        reactants_3a = Nucleus.cast_list(["he4", "he4", "he4"])
        products_3a = [Nucleus("c12")]

        # Get all available RatePairs, i.e. rates that have reverse rates.
        # It also filters out removed rates since it considers self.rates
        # Then do additional filtering based on different criteria
        for rp in self.get_rate_pairs():
            reactants = rp.forward.reactants
            products = rp.forward.products
            if rp.forward.weak or rp.reverse.weak:
                continue
            if len(reactants) > 2 or len(products) > 2:
                if reactants == reactants_3a and products == products_3a:
                    nse_rate_pairs.append(rp)
                continue

            non_LIG_count = sum(nuc not in LIG for nuc in reactants + products)
            if non_LIG_count in (1, 2):
                nse_rate_pairs.append(rp)

        return nse_rate_pairs

    def _fill_npa_index(self, n_indent, of):
        """Write 1-based indices of light isotopes (p, n, he4) if present in the network
        (set to -1 if absent).

        """

        LIG = Nucleus.cast_list(["p", "n", "he4"])
        for nuc in LIG:
            idx = self.unique_nuclei.index(nuc) + 1 if nuc in self.unique_nuclei else -1
            of.write(f"{self.indent*n_indent}constexpr int {nuc.short_spec_name.capitalize()}_index = {idx};\n")

    def _write_nse_rate_pair_data(self, n_indent, of):
        """Write the RatePair data table, 2D array of shape (NumNSERatePairs, 8),
        needed for the NSE_NET algorithm.

        Each row corresponds to a RatePair returned by _get_nse_rate_pairs().
        Nuclei indices follow NetworkSpecies and rate indices follow NetworkRates.

        - Columns 1–3: 1-based index of reactant nuclei (forward rate)
        - Columns 4–6: 1-based index of product nuclei (forward rate)
        - Column 7:    1-based index of the forward rate
        - Column 8:    1-based index of the reverse rate

        """

        nse_rate_pairs = self._get_nse_rate_pairs()
        NumNSERatePairs = len(nse_rate_pairs)
        dtype = _signed_rate_dtype(NumNSERatePairs)

        # Write Fill in the rate indices
        of.write(f"{self.indent*n_indent}constexpr int NumNSERatePairs = {NumNSERatePairs};\n\n")
        of.write(f"{self.indent*n_indent}inline AMREX_GPU_MANAGED amrex::Array2D<{dtype}, 1, NumNSERatePairs, 1, 8, amrex::Order::C> rate_pair_data {{\n")

        for n, rp in enumerate(nse_rate_pairs):
            fr = rp.forward
            rr = rp.reverse

            # Find the reactants and products indices for the forward rate
            # Change nuclei indices to 1-based
            reactant_idx = [-1 for n in range(3 - len(fr.reactants))]
            product_idx = [-1 for n in range(3 - len(fr.products))]

            for nuc in fr.reactants:
                reactant_idx.append(self.unique_nuclei.index(nuc) + 1)
            for nuc in fr.products:
                product_idx.append(self.unique_nuclei.index(nuc) + 1)

            reactant_idx.sort()
            product_idx.sort()

            # Find rate index and note that they are 1-based
            fr_idx = self.all_rates.index(fr) + 1
            rr_idx = self.all_rates.index(rr) + 1

            of.write(f"{self.indent*(n_indent+1)}"
                     f"{reactant_idx[0]}, {reactant_idx[1]}, {reactant_idx[2]}, "
                     f"{product_idx[0]}, {product_idx[1]}, {product_idx[2]}, "
                     f"{fr_idx}, {rr_idx}")

            if n < NumNSERatePairs - 1:
                of.write(",")
            of.write(f"  // fr: {fr.fname}, rr: {rr.fname}\n")

        of.write(f"{self.indent*n_indent}}};\n")
