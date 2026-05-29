"""A simple C++ reaction network for integrating into other C++ codes"""


from pathlib import Path

from pynucastro.networks.base_cxx_network import BaseCxxNetwork
from pynucastro.screening import get_screening_pair_set


class SimpleCxxNetwork(BaseCxxNetwork):
    """A basic C++ network."""

    def __init__(self, *args, **kwargs):

        # Initialize BaseCxxNetwork parent class
        super().__init__(*args, **kwargs)

        self.function_specifier = "inline"
        self.gpu_data_specifier = ""
        self.dtype = "Real"
        self.array_namespace = ""

    def _get_template_files(self):

        path = self.pynucastro_dir/"templates/simple-cxx-network"

        return path.glob("*.template")

    def _compute_all_screening_factors(self, n_indent, of):
        """Composes the screening factors string for all rates.
        It evaluates log(screening) and stores them to rate_eval.log_screen.
        It doesn't include the Temperature derivative since we currently
        do not support temperature evolution in simple-cxx network.

        """
        self._compute_screening_factors(self, n_indent, of, rate=self.get_rates(),
                                        do_T_derivatives=False)

    def _write_network(self, odir=None):
        """Output the RHS, jacobian and ancillary files for the
        system of ODEs that this network describes, using the template
        files.

        """

        assert len(self.temperature_tabular_rates) == 0, "SimpleCxxNetwork does not support TemperatureTabular rates"
        assert len(self.starlib_rates) == 0, "SimpleCxxNetwork does not support StarLib rates"

        super()._write_network(odir=odir)

        if odir is None:
            odir = Path.cwd()
        # create a header file with the nuclei properties
        with open(Path(odir, "network_properties.H"), "w") as of:
            of.write("#ifndef NETWORK_PROPERTIES_H\n")
            of.write("#define NETWORK_PROPERTIES_H\n")
            of.write("#include <vector>\n")
            of.write("#include <string>\n")
            of.write("#include <amrex_bridge.H>\n\n")

            of.write(f"constexpr int NumSpec = {len(self.unique_nuclei)};\n\n")
            of.write(f"constexpr int NumSpecExtra = {len(self.approx_nuclei)};\n\n")
            of.write("constexpr int NumSpecTotal = NumSpec + NumSpecExtra;\n\n")

            of.write("// Note: these are 0-based\n")

            of.write("constexpr Real aion[NumSpecTotal] = {\n")
            for n, nuc in enumerate(self.unique_nuclei + self.approx_nuclei):
                of.write(f"    {nuc.A:6.1f}, // {n} : {nuc}\n")
            of.write(" };\n\n")

            of.write("constexpr Real aion_inv[NumSpecTotal] = {\n")
            for n, nuc in enumerate(self.unique_nuclei + self.approx_nuclei):
                of.write(f"    1.0/{nuc.A:6.1f}, // {n} : {nuc}\n")
            of.write(" };\n\n")

            of.write("constexpr Real zion[NumSpecTotal] = {\n")
            for n, nuc in enumerate(self.unique_nuclei + self.approx_nuclei):
                of.write(f"    {nuc.Z:6.1f}, // {n} : {nuc}\n")
            of.write(" };\n\n")

            of.write("static const std::vector<std::string> spec_names = {\n")
            for n, nuc in enumerate(self.unique_nuclei + self.approx_nuclei):
                of.write(f"    \"{nuc.short_spec_name.capitalize()}\", // {n}\n")
            of.write(" };\n\n")

            of.write("namespace Species {\n")
            of.write("  enum NetworkSpecies {\n")
            for n, nuc in enumerate(self.unique_nuclei + self.approx_nuclei):
                if n == 0:
                    of.write(f"    {nuc.short_spec_name.capitalize()}=1,\n")
                else:
                    of.write(f"    {nuc.short_spec_name.capitalize()},\n")
            of.write("  };\n")
            of.write("}\n\n")

            of.write("#endif\n")
