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
            for nuc in self.unique_nuclei:
                of.write(f"{nuc.spec_name:25} {nuc.short_spec_name:6} {nuc.A:6.1f} {nuc.Z:6.1f}\n")

            for nuc in self.approx_nuclei:
                of.write(f"__extra_{nuc.spec_name:17} {nuc.short_spec_name:6} {nuc.A:6.1f} {nuc.Z:6.1f}\n")
