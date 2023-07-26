"""A simple C++ reaction network for integrating into other C++ codes"""


import glob
import os

from pynucastro.networks.base_cxx_network import BaseCxxNetwork
from pynucastro.nucdata import Nucleus
from pynucastro.rates import ReacLibRate


class SimpleCxxNetwork(BaseCxxNetwork):
    def __init__(self, *args, **kwargs):

        # Initialize BaseCxxNetwork parent class
        super().__init__(*args, **kwargs)

        self.function_specifier = "inline"
        self.dtype = "Real"

    def _get_template_files(self):

        template_pattern = os.path.join(self.pynucastro_dir,
                                        'templates',
                                        'simple-cxx-network',
                                        '*.template')

        return glob.glob(template_pattern)

    def _write_network(self, odir=None):
        """
        This writes the RHS, jacobian and ancillary files for the system of ODEs that
        this network describes, using the template files.
        """

        super()._write_network(odir=odir)

        if odir is None:
            odir = os.getcwd()
        # create a header file with the nuclei properties
        with open(os.path.join(odir, "network_properties.H"), "w") as of:
            of.write("#ifndef NETWORK_PROPERTIES_H\n")
            of.write("#define NETWORK_PROPERTIES_H\n")
            of.write("#include <vector>\n")
            of.write("#include <string>\n")
            of.write("#include <amrex_bridge.H>\n\n")

            of.write(f"constexpr int NumSpec = {len(self.unique_nuclei)};\n\n")

            of.write("constexpr Real aion[NumSpec] = {\n")
            for n, nuc in enumerate(self.unique_nuclei):
                of.write(f"    {nuc.A:6.1f}, // {n}\n")
            of.write(" };\n\n")

            of.write("constexpr Real aion_inv[NumSpec] = {\n")
            for n, nuc in enumerate(self.unique_nuclei):
                of.write(f"    1.0/{nuc.A:6.1f}, // {n}\n")
            of.write(" };\n\n")

            of.write("constexpr Real zion[NumSpec] = {\n")
            for n, nuc in enumerate(self.unique_nuclei):
                of.write(f"    {nuc.Z:6.1f}, // {n}\n")
            of.write(" };\n\n")

            of.write("static const std::vector<std::string> spec_names = {\n")
            for n, nuc in enumerate(self.unique_nuclei):
                of.write(f"    \"{nuc.short_spec_name.capitalize()}\", // {n}\n")
            of.write(" };\n\n")

            of.write("namespace Species {\n")
            of.write("  enum NetworkSpecies {\n")
            for n, nuc in enumerate(self.unique_nuclei):
                if n == 0:
                    of.write(f"    {nuc.short_spec_name.capitalize()}=1,\n")
                else:
                    of.write(f"    {nuc.short_spec_name.capitalize()},\n")
            of.write("  };\n")
            of.write("}\n\n")

            of.write("#endif\n")
