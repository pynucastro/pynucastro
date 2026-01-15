"""Classes and methods for managing a tabular electron / positron EOS.
This uses a tabulation of the Helmholtz free energy and interpolates
on it to find the necessary thermodynamic quantities.
"""

import bz2
import re

import numpy as np

from pynucastro.constants import constants
from pynucastro.rates.files import _pynucastro_dir


class Helmholtz:
    """a container to hold the Helmholtz state and its derivatives"""

    __slots__ = ("F", "dF_drho", "dF_dT",
                 "d2F_drho2", "d2F_dT2", "d2F_drhodT",
                 "d3F_drho2dT", "d3F_drhodT2", "d4F_drho2dT2")

    def __init__(self,
                 F=None, dF_drho=None, dF_dT=None,
                 d2F_drho2=None, d2F_dT2=None, d2F_drhodT=None,
                 d3F_drho2dT=None, d3F_drhodT2=None, d4F_drho2dT2=None):
        self.F = float(F)
        self.dF_drho = float(dF_drho)
        self.dF_dT = float(dF_dT)
        self.d2F_drho2 = float(d2F_drho2)
        self.d2F_dT2 = float(d2F_dT2)
        self.d2F_drhodT = float(d2F_drhodT)
        self.d3F_drho2dT = float(d3F_drho2dT)
        self.d3F_drhodT2 = float(d3F_drhodT2)
        self.d4F_drho2dT2 = float(d4F_drho2dT2)


class ThermoQuantity:
    """a container to hold a thermodynamic quantity, its first
    derivatives, and its 2nd crossed derivative."""

    __slots__ = ("q", "dq_drho", "dq_dT", "d2q_drhodT")

    def __init__(self, q=None, dq_drho=None, dq_dT=None, d2q_drhodT=None):
        self.q = float(q)
        self.dq_drho = float(dq_drho)
        self.dq_dT = float(dq_dT)
        self.d2q_drhodT = float(d2q_drhodT)


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
                        self.log10_rho_min = float(g.group(1))
                        self.log10_rho_max = float(g.group(2))
                        self.rho_npts = int(g.group(3))
                    elif g := temp_line.match(line):
                        self.log10_temp_min = float(g.group(1))
                        self.log10_temp_max = float(g.group(2))
                        self.temp_npts = int(g.group(3))
                    elif g := version_line.match(line):
                        self.table_version = g.group(1)

        # now read and store the data in 4 separate lists
        with bz2.open(_pynucastro_dir / "eos" / self.table, mode="rt") as tf:

            # the first rho_npts * temp_npts lines are the Helmholtz data
            self.helm = []
            for _ in range(self.rho_npts * self.temp_npts):
                line = tf.readline()
                fields = line.split()
                h = Helmholtz(F=fields[0],
                              dF_drho=fields[1], dF_dT=fields[2],
                              d2F_drho2=fields[3], d2F_dT2=fields[4], d2F_drhodT=fields[5],
                              d3F_drho2dT=fields[6], d3F_drhodT2=fields[7], d4F_drho2dT2=fields[8])
                self.helm.append(h)

            # the next rho_npts * temp_npts lines are pressure derivatives
            self.dp_drho = []
            for _ in range(self.rho_npts * self.temp_npts):
                line = tf.readline()
                fields = line.split()
                q = ThermoQuantity(q=fields[0],
                                   dq_drho=fields[1], dq_dT=fields[2],
                                   d2q_drhodT=fields[3])
                self.dp_drho.append(q)

            # the next rho_npts * temp_npts lines are the degeneracy parameter
            self.psi = []
            for _ in range(self.rho_npts * self.temp_npts):
                line = tf.readline()
                fields = line.split()
                q = ThermoQuantity(q=fields[0],
                                   dq_drho=fields[1], dq_dT=fields[2],
                                   d2q_drhodT=fields[3])
                self.psi.append(q)

            # the final rho_npts * temp_npts lines are the number density
            self.ne = []
            for _ in range(self.rho_npts * self.temp_npts):
                line = tf.readline()
                fields = line.split()
                q = ThermoQuantity(q=fields[0],
                                   dq_drho=fields[1], dq_dT=fields[2],
                                   d2q_drhodT=fields[3])
                self.ne.append(q)

    def _index(self, irho, jtemp):
        # in our 1D tabulation, density varies the fastest
        return jtemp * self.rho_npts + irho

    def __str__(self):
        ostr = ""
        ostr += "electron/positron tabular EOS\n"
        ostr += f"temperature range: [10**{self.log10_temp_min}, 10**{self.log10_temp_max}] with {self.temp_npts} points\n"
        ostr += f"density range: [10**{self.log10_rho_min}, 10**{self.log10_rho_max}] with {self.rho_npts} points\n"
        ostr += f"table generation git version: {self.table_version}\n"
        return ostr


