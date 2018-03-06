"""A Fortran reaction network for integration into the StarKiller
Microphysics set of reaction networks used by astrophysical hydrodynamics
codes"""

from __future__ import print_function

import glob
import os

from pynucastro.networks import BaseFortranNetwork

class StarKillerNetwork(BaseFortranNetwork):
    def __init__(self, *args, **kwargs):
        # Initialize BaseFortranNetwork parent class
        super(StarKillerNetwork, self).__init__(*args, **kwargs)

        # Set up some directories
        self.template_dir = os.path.join(self.pynucastro_dir,
                                         'templates',
                                         'starkiller-microphysics')
        self.template_file_select = os.path.join(self.template_dir,
                                                 '*.template')
        self.template_files = glob.glob(self.template_file_select)

    def _initial_mass_fractions(self, n_indent, of):
        # Redefine initial mass fractions tag to set the
        # mass fractions in the burn_cell unit test inputs file.
        for i, n in enumerate(self.unique_nuclei):
            of.write("\n# {nuc: <5} initial mass fraction\n".format(nuc=str(n)))
            of.write("{}massfractions({}) = 0.0d0\n".format(
                self.indent*n_indent, i))
