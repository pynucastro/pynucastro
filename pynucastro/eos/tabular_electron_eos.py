"""Classes and methods for managing a tabular electron / positron EOS.
This uses a tabulation of the Helmholtz free energy and interpolates
on it to find the necessary thermodynamic quantities.
"""

import bz2
import re

import numpy as np

from pynucastro.constants import constants
from pynucastro.rates.files import _pynucastro_dir


class TabularElectronEOS:
    """tabular EOS"""

    def __init__(self, table="helm_table_p128_q800.dat.bz2"):
        self.table = table

        # read the table in two passes.  First the metadata at the end.
        rho_line = re.compile(r"# log10\(rho_lo\) = ([-\.\d]*), log10\(rho_hi\) = ([-\.\d]*), rho_pts = ([\d]*)")
        temp_line = re.compile(r"# log10\(T_lo\) = ([-\.\d]*), log10\(T_hi\) = ([-\.\d]*), T_pts = ([\d]*)")
        version_line = re.compile(r"# git version: ([\w\.\-]+)")

        with bz2.open(_pynucastro_dir / "eos" / self.table, mode="rt") as tf:
            while line := tf.readline():
                if line.startswith("#"):
                    if g := rho_line.match(line):
                        self.log10_rho_min = g.group(1)
                        self.log10_rho_max = g.group(2)
                        self.rho_npts = g.group(3)
                    elif g := temp_line.match(line):
                        self.log10_temp_min = g.group(1)
                        self.log10_temp_max = g.group(2)
                        self.temp_npts = g.group(3)
                    elif g := version_line.match(line):
                        self.table_version = g.group(1)

    def __str__(self):
        ostr = ""
        ostr += "electron/positron tabular EOS\n"
        ostr += f"temperature range: [10**{self.log10_temp_min}, 10**{self.log10_temp_max}] with {self.temp_npts} points\n"
        ostr += f"density range: [10**{self.log10_rho_min}, 10**{self.log10_rho_max}] with {self.rho_npts} points\n"
        ostr += f"table generation git version: {self.table_version}\n"
        return ostr


