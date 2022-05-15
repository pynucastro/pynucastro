"""A Fortran reaction network for integration into the StarKiller
Microphysics set of reaction networks used by astrophysical hydrodynamics
codes"""


import glob
import os

from pynucastro.networks import BaseCxxNetwork


class StarKillerCxxNetwork(BaseCxxNetwork):
    def __init__(self, rate_params=None, *args, **kwargs):
        # Initialize BaseFortranNetwork parent class
        super().__init__(*args, **kwargs)

        self.rate_params = rate_params

    def _get_template_files(self):

        template_pattern = os.path.join(self.pynucastro_dir,
                                        'templates',
                                        'starkiller-cxx-microphysics',
                                        '*.template')

        return glob.glob(template_pattern)

    def _write_network(self, odir=None):
        """
        This writes the RHS, jacobian and ancillary files for the system of ODEs that
        this network describes, using the template files.
        """

        super()._write_network(odir=odir)

        # create a .net file with the nuclei properties
        with open("pynucastro.net", "w") as of:
            for nuc in self.unique_nuclei:
                of.write("{:25} {:6} {:6.1f} {:6.1f}\n".format(
                    nuc.spec_name, nuc.short_spec_name, nuc.A, nuc.Z))

        # write out some network properties
        with open("NETWORK_PROPERTIES", "w") as of:
            of.write(f"NSCREEN := {self.num_screen_calls}\n")

        # write the _parameters file
        with open("_parameters", "w") as of:
            of.write("@namespace: network\n\n")
            if self.rate_params:
                for p in self.rate_params:
                    of.write(f"{p}    int     0\n")
