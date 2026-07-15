"""Classes and methods for managing a tabular electron / positron EOS.
This uses a tabulation of the Helmholtz free energy and interpolates
on it to find the necessary thermodynamic quantities.
"""

import bz2
import re

import numpy as np

from pynucastro.rates.files import _pynucastro_dir
from .eos_components import EOSComponentState


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


def psi0(z):
    """The ψ0 Hermite basis function and its derivatives, TS00 Eq. 15"""
    _psi0 = z**3 * (z * (-6.0 * z + 15.0) - 10.0) + 1.0
    _dpsi0_dz = z**2 * (z * (-30.0 * z + 60.0) - 30.0)
    _d2psi0_dz2 = z * (z * (-120.0e0 * z + 180.0) - 60.0)
    return _psi0, _dpsi0_dz, _d2psi0_dz2


def psi1(z):
    """The ψ1 Hermite basis function and its derivatives, TS00 Eq. 15"""
    _psi1 = z * (z**2 * (z * (-3.0 * z + 8.0) - 6.0) + 1.0)
    _dpsi1_dz = z**2 * (z * (-15.0 * z + 32.0) - 18.0) + 1.0
    _d2psi1_dz2 = z * (z * (-60.0 * z + 96.0) - 36.0)
    return _psi1, _dpsi1_dz, _d2psi1_dz2


def psi2(z):
    """The ψ2 Hermite basis function and its derivatives, TS00 Eq. 15"""
    _psi2 = 0.5 * z**2 * (z * (z * (-z + 3.0) - 3.0) + 1.0)
    _dpsi2_dz = 0.5 * z * (z * (z * (-5.0 * z + 12.0) - 9.0) + 2.0)
    _d2psi2_dz2 = 0.5 * (z * (z * (-20.0 * z + 36.0) - 18.0) + 2.0)
    return _psi2, _dpsi2_dz, _d2psi2_dz2


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

        # store spacings
        self.dlogrho = (self.log10_rho_max - self.log10_rho_min) / float(self.rho_npts - 1)
        self.dlogrho_inv = 1.0 / self.dlogrho

        self.dlogT = (self.log10_temp_max - self.log10_temp_min) / float(self.temp_npts - 1)
        self.dlogT_inv = 1.0 / self.dlogT

        # compute the density and temperature grid points
        self.rho = np.logspace(self.log10_rho_min, self.log10_rho_max, self.rho_npts, endpoint=True)
        self.T = np.logspace(self.log10_temp_min, self.log10_temp_max, self.temp_npts, endpoint=True)

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

    def interp_helmholtz(self, i, j):
        pass

    def pe_state(self, rho=None, T=None, comp=None):
        """Find the pressure and energy given density, temperature,
        and composition"""

        # find the i and j location (density and temperature indices into the table)
        rho_e = rho * comp.ye()
        _iat = int((np.log10(rho_e) - self.log10_rho_min) * self.dlogrho_inv) + 1
        iat = np.clip(_iat, 1, self.rho_npts-1) - 1

        _jat = int((np.log10(T) - self.log10_temp_min) * self.dlogT_inv) + 1
        jat = np.clip(_jat, 1, self.temp_npts-1) - 1

        # compute the dimensionless quantities used in the basis functions, TS00 Eq. 18
        drho = self.rho[iat+1] - self.rho[iat]
        dT = self.T[jat+1] - self.T[jat]

        x = max((rho_e - self.rho[iat]) / drho, 0.0)
        y = max((T - self.T[jat]) / dT, 0.0)

        # get the basis functions
        psi0_x, dpsi0_x, d2psi0_x = psi0(x)
        psi0_1mx, dpsi0_1mx, d2psi0_1mx = psi0(1.0 - x)

        psi1_x, dpsi1_x, d2psi1_x = psi1(x)
        psi1_1mx, dpsi1_1mx, d2psi1_1mx = psi1(1.0 - x)

        psi2_x, dpsi2_x, d2psi2_x = psi2(x)
        psi2_1mx, dpsi2_1mx, d2psi2_1mx = psi2(1.0 - x)

        psi0_y, dpsi0_y, d2psi0_y = psi0(y)
        psi0_1my, dpsi0_1my, d2psi0_1my = psi0(1.0 - y)

        psi1_y, dpsi1_y, d2psi1_y = psi1(y)
        psi1_1my, dpsi1_1my, d2psi1_1my = psi1(1.0 - y)

        psi2_y, dpsi2_y, d2psi2_y = psi2(y)
        psi2_1my, dpsi2_1my, d2psi2_1my = psi2(1.0 - y)

        
        
        ele_state = EOSComponentState(eta=eta,
                                      n=n_e, p=p_e, e=e_e,
                                      dn_drho=dne_drho, dn_dT=dne_dT,
                                      dp_drho=dpe_drho, dp_dT=dpe_dT,
                                      de_drho=dee_drho, de_dT=dee_dT)



    def __str__(self):
        ostr = ""
        ostr += "electron/positron tabular EOS\n"
        ostr += f"temperature range: [10**{self.log10_temp_min}, 10**{self.log10_temp_max}] with {self.temp_npts} points\n"
        ostr += f"density range: [10**{self.log10_rho_min}, 10**{self.log10_rho_max}] with {self.rho_npts} points\n"
        ostr += f"table generation git version: {self.table_version}\n"
        return ostr


