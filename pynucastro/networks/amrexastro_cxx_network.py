"""A C++ reaction network for integration into the AMReX Astro
Microphysics set of reaction networks used by astrophysical hydrodynamics
codes"""


import glob
import os

from pynucastro.networks.base_cxx_network import BaseCxxNetwork


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
                of.write(f"{self.indent*n_indent}if (disable_{r.fname}) {{\n")
                of.write(f"{self.indent*n_indent}    rate_eval.screened_rates(k_{r.fname}) = 0.0;\n")
                of.write(f"{self.indent*n_indent}    if constexpr (std::is_same<T, rate_derivs_t>::value) {{\n")
                of.write(f"{self.indent*n_indent}        rate_eval.dscreened_rates_dT(k_{r.fname}) = 0.0;\n")
                of.write(f"{self.indent*n_indent}    }}\n")
                # check for the reverse too -- we disable it with the same parameter
                rr = self.find_reverse(r)
                if rr is not None:
                    of.write(f"{self.indent*n_indent}    rate_eval.screened_rates(k_{rr.fname}) = 0.0;\n")
                    of.write(f"{self.indent*n_indent}    if constexpr (std::is_same<T, rate_derivs_t>::value) {{\n")
                    of.write(f"{self.indent*n_indent}    rate_eval.dscreened_rates_dT(k_{rr.fname}) = 0.0;\n")
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

        # write out some network properties
        with open(os.path.join(odir, "NETWORK_PROPERTIES"), "w") as of:
            of.write(f"NSCREEN := {self.num_screen_calls}\n")

        # write the _parameters file
        with open(os.path.join(odir, "_parameters"), "w") as of:
            of.write("@namespace: network\n\n")
            if self.disable_rate_params:
                for r in self.disable_rate_params:
                    of.write(f"disable_{r.fname}    int     0\n")
