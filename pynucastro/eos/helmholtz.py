import numpy as np
from interface import EosState

""" Here we implement the helmholtz eos. The idea is to call for two thermodynamic variables and the composition
of a reaction network."""

class HelmholtzTable:

    #The first step is to interpolate the electron_positron eos table by the use of biquintic
    #polinomials.

    def __init__(self, table='helm_table.dat'):

        self.table=table

        self.jmax = 201   #201 number of temp pts
        self.tlo = 3.0
        self.thi = 13.0
        self.tstp  = (self.thi - self.tlo) / (self.jmax-1)
        self.tstpi = 1.0 / self.tstp

        self.imax = 541   #541 number of dens pts
        self.dlo   = -12.0
        self.dhi   = 15.0
        self.dstp  = (self.dhi - self.dlo) / (self.imax-1)
        self.dstpi = 1.0 / self.dstp

        self.f = None
        self.pd = None
        self.ef = None
        self.xnf = None

        self._read_file()

        self.din = None
        self.T = None
        self.ye = None

        self.iat = None
        self.jat = None

    # Now, let us define a set of indices on which the table will
    # operate
    #=============================================================
    def _rho_T_indices(self):

        iat = int((np.log10(self.din) - self.dlo) * self.dstpi)
        iat = max(0, min(iat, self.imax-2))

        jat = int((np.log10(self.T) - self.tlo) * self.tstpi)
        jat = max(0, min(jat, self.jmax-2))

        return iat, jat

    def _rho_from_index(self, i):

        logrho_i = self.dlo + i * self.dstp
        return 10.0**logrho_i

    def _T_from_index(self, j):

        logT_j = self.tlo + j * self.tstp
        return 10.0**logT_j

    def _read_file(self):
        """ From here we read the table only once to interpolate F"""

        def _extract_line_data(line):
            line_entries = line.strip().split()
            f = np.array([np.float64(s) for s in line_entries])
            return f

        with open(self.table, 'r') as file:

            f_local = np.empty((self.imax, self.jmax, 9))
            for j in range(0,self.jmax):
                for i in range(0,self.imax):
                    line = file.readline()
                    f_local[i,j,:] = _extract_line_data(line)[:]

            pd_local = np.empty((self.imax, self.jmax, 4))
            for j in range(0,self.jmax):
                for i in range(0,self.imax):
                    line = file.readline()
                    pd_local[i,j,:] = _extract_line_data(line)[:]

            ef_local = np.empty((self.imax, self.jmax, 4))
            for j in range(0,self.jmax):
                for i in range(0,self.imax):
                    line = file.readline()
                    ef_local[i,j,:] = _extract_line_data(line)[:]

            xnf_local = np.empty((self.imax, self.jmax, 4))
            for j in range(0,self.jmax):
                for i in range(0,self.imax):
                    line = file.readline()
                    xnf_local[i,j, :] = _extract_line_data(line)[:]

        self.f = f_local
        self.pd = pd_local
        self.ef = ef_local
        self.xnf = xnf_local

    #Let us define the bicubic polynomials:
    #======================================

    # xpsi0
    def _xpsi0(self, z):

        return z * z * (2.0 * z - 3.0) + 1.0

    def _xdpsi0(self, z):

        return z * (6.0 * z - 6.0)

    # xpsi1
    def _xpsi1(self, z):

        return z * (z * (z - 2.0) + 1.0)

    def _xdpsi1(self, z):

        return z * (3.0 * z - 4.0) + 1.0

    #Let us define the biquintic the polynomials:
    #============================================

    # psi0
    def _psi0(self, z):

        return  z * z * z * (z * (-6.0 * z + 15.0) -10.0) + 1.0

    def _dpsi0(self,z):

        return z * z * (z * (-30.0*z + 60.0) - 30.0)

    def _ddpsi0(self,z):

        return z * (z * (-120.0 * z + 180.0) - 60.0)

    # psi1
    def _psi1(self, z):

        return z * ( z * z * (z * (-3.0 * z + 8.0) - 6.0) + 1.0)

    def _dpsi1(self,z):

        return z * z * ( z * (-15.0 * z + 32.0) - 18.0) + 1.0

    def _ddpsi1(self,z):

        return z * (z * (-60.0 * z + 96.0) - 36.0)

    # psi2
    def _psi2(self, z):

        return 0.5 * z * z * ( z * ( z * (-z + 3.0) - 3.0) + 1.0)

    def _dpsi2(self,z):

        return 0.5 * z * (z * (z * (-5.0 * z + 12.0) - 9.0) + 2.0)

    def _ddpsi2(self,z):

        return 0.5 * (z * ( z * (-20.0 * z + 36.0) - 18.0) + 2.0)

    # Let us read the density coeficients that participate in the
    # interpolation:
    #============================================================

    def _read_wd(self, der=0, sel='biquintic'):

        x = (self.din - self._rho_from_index(self.iat)) / \
            (self._rho_from_index(self.iat+1) - self._rho_from_index(self.iat))

        def _no_derivative_coeff(self, x):

            psi = np.zeros(6)
            din_coeff = np.zeros(6)
            wd = np.zeros(6)

            d_din = self._rho_from_index(self.iat+1) - self._rho_from_index(self.iat)
            d2_din = (self._rho_from_index(self.iat+1) - self._rho_from_index(self.iat))**2

            psi[0] = self._psi0(x)
            psi[1] = self._psi0(1-x)
            psi[2] = self._psi1(x)
            psi[3] = self._psi1(1-x)
            psi[4] = self._psi2(x)
            psi[5] = self._psi2(1-x)

            din_coeff[0] = 1.0
            din_coeff[1] = 1.0
            din_coeff[2] = d_din
            din_coeff[3] = -d_din
            din_coeff[4] = d2_din
            din_coeff[5] = d2_din

            wd = psi * din_coeff
            return wd

        def _first_derivative_coeff(self, x):

            psi = np.zeros(6)
            din_coeff = np.zeros(6)
            wd = np.zeros(6)

            d_din = self._rho_from_index(self.iat+1) - self._rho_from_index(self.iat)
            d_din_i = 1.0 / d_din

            psi[0] = self._dpsi0(x)
            psi[1] = self._dpsi0(1-x)
            psi[2] = self._dpsi1(x)
            psi[3] = self._dpsi1(1-x)
            psi[4] = self._dpsi2(x)
            psi[5] = self._dpsi2(1-x)

            din_coeff[0] = d_din_i
            din_coeff[1] = -d_din_i
            din_coeff[2] = 1.0
            din_coeff[3] = 1.0
            din_coeff[4] = d_din
            din_coeff[5] = -d_din

            wd = psi * din_coeff
            return wd

        def _second_derivative_coeff(self, x):

            psi = np.zeros(6)
            din_coeff = np.zeros(6)
            wd = np.zeros(6)

            d_din = self._rho_from_index(self.iat+1) - self._rho_from_index(self.iat)
            d_din_i = 1.0 / d_din

            d2_din = (self._rho_from_index(self.iat+1) - self._rho_from_index(self.iat))**2
            d2_din_i = 1.0 / d2_din

            psi[0] = self._ddpsi0(x)
            psi[1] = self._ddpsi0(1-x)
            psi[2] = self._ddpsi1(x)
            psi[3] = self._ddpsi1(1-x)
            psi[4] = self._ddpsi2(x)
            psi[5] = self._ddpsi2(1-x)

            din_coeff[0] = d2_din_i
            din_coeff[1] = d2_din_i
            din_coeff[2] = d_din_i
            din_coeff[3] = -d_din_i
            din_coeff[4] = 1.0
            din_coeff[5] = 1.0

            wd = psi * din_coeff
            return wd

        def _no_xderivative_coeff(self, x):

            psi = np.zeros(4)
            din_coeff = np.zeros(4)
            wd = np.zeros(4)

            d_din = self._rho_from_index(self.iat+1) - self._rho_from_index(self.iat)

            psi[0] = self._psi0(x)
            psi[1] = self._psi0(1-x)
            psi[2] = self._psi1(x)
            psi[3] = self._psi1(1-x)

            din_coeff[0] = 1.0
            din_coeff[1] = 1.0
            din_coeff[2] = d_din
            din_coeff[3] = -d_din

            wd = psi * din_coeff
            return wd

        def _first_xderivative_coeff(self, x):

            psi = np.zeros(4)
            din_coeff = np.zeros(4)
            wd = np.zeros(4)

            d_din = self._rho_from_index(self.iat+1) - self._rho_from_index(self.iat)
            d_din_i = 1.0 / d_din

            psi[0] = self._xdpsi0(x)
            psi[1] = self._xdpsi0(1-x)
            psi[2] = self._xdpsi1(x)
            psi[3] = self._xdpsi1(1-x)

            din_coeff[0] = d_din_i
            din_coeff[1] = -d_din_i
            din_coeff[2] = 1.0
            din_coeff[3] = 1.0

            wd = psi * din_coeff
            return wd

        if sel == "biquintic":

            wd = np.zeros(6)
            coeff_derivatives = {0: _no_derivative_coeff,
                                 1: _first_derivative_coeff,
                                 2: _second_derivative_coeff}

        elif sel == 'bicubic':

            wd = np.zeros(4)
            coeff_derivatives = {0: _no_xderivative_coeff,
                                 1: _first_xderivative_coeff}
        else:

            raise NotImplementedError("There is no implementation for the selected type")

        func = coeff_derivatives.get(der)
        wd = func(self, x)

        return wd

    # Now we need to construct the temperature coefficients that participate
    # in the table interpolation.

    def _read_wt(self, der=0, sel='biquintic'):

        y = (self.T - self._T_from_index(self.jat)) / \
            (self._T_from_index(self.jat+1) - self._T_from_index(self.jat))

        def _no_derivative_coeff(self, y):

            tin_coeff = np.zeros(6)
            psi = np.zeros(6)
            wt = np.zeros(6)

            d_tin = self._T_from_index(self.jat+1) - self._T_from_index(self.jat)
            d2_tin = (self._T_from_index(self.jat+1) - self._T_from_index(self.jat))**2

            psi[0] = self._psi0(y)
            psi[1] = self._psi0(1-y)
            psi[2] = self._psi1(y)
            psi[3] = self._psi1(1-y)
            psi[4] = self._psi2(y)
            psi[5] = self._psi2(1-y)

            tin_coeff[0] = 1.0
            tin_coeff[1] = 1.0
            tin_coeff[2] = d_tin
            tin_coeff[3] = -d_tin
            tin_coeff[4] = d2_tin
            tin_coeff[5] = d2_tin

            wt = psi * tin_coeff
            return wt

        def _first_derivative_coeff(self, y):

            tin_coeff = np.zeros(6)
            psi = np.zeros(6)
            wt = np.zeros(6)

            d_tin = self._T_from_index(self.jat+1) - self._T_from_index(self.jat)
            d_tin_i = 1.0 / d_tin

            psi[0] = self._dpsi0(y)
            psi[1] = self._dpsi0(1-y)
            psi[2] = self._dpsi1(y)
            psi[3] = self._dpsi1(1-y)
            psi[4] = self._dpsi2(y)
            psi[5] = self._dpsi2(1-y)

            tin_coeff[0] = d_tin_i
            tin_coeff[1] = -d_tin_i
            tin_coeff[2] = 1.0
            tin_coeff[3] = 1.0
            tin_coeff[4] = d_tin
            tin_coeff[5] = -d_tin

            wt = psi * tin_coeff
            return wt

        def _second_derivative_coeff(self, y):

            tin_coeff = np.zeros(6)
            psi = np.zeros(6)
            wt = np.zeros(6)

            d_tin = self._T_from_index(self.jat+1) - self._T_from_index(self.jat)
            d_tin_i = 1.0 / d_tin

            d2_tin = (self._T_from_index(self.jat+1) - self._T_from_index(self.jat))**2
            d2_tin_i = 1.0 / d2_tin

            psi[0] = self._ddpsi0(y)
            psi[1] = self._ddpsi0(1-y)
            psi[2] = self._ddpsi1(y)
            psi[3] = self._ddpsi1(1-y)
            psi[4] = self._ddpsi2(y)
            psi[5] = self._ddpsi2(y)

            tin_coeff[0] = d2_tin_i
            tin_coeff[1] = d2_tin_i
            tin_coeff[2] = d_tin_i
            tin_coeff[3] = -d_tin_i
            tin_coeff[4] = 1.0
            tin_coeff[5] = 1.0

            wt = psi * tin_coeff
            return wt

        def _no_xderivative_coeff(self, y):

            psi = np.zeros(4)
            din_coeff = np.zeros(4)
            wt = np.zeros(4)

            d_tin = self._T_from_index(self.jat+1) - self._T_from_index(self.jat)

            psi[0] = self._psi0(y)
            psi[1] = self._psi0(1-y)
            psi[2] = self._psi1(y)
            psi[3] = self._psi1(1-y)

            din_coeff[0] = 1.0
            din_coeff[1] = 1.0
            din_coeff[2] = d_tin
            din_coeff[3] = -d_tin

            wt = psi * din_coeff
            return wt

        def _first_xderivative_coeff(self, y):

            psi = np.zeros(4)
            din_coeff = np.zeros(4)
            wt = np.zeros(4)

            d_tin = self._T_from_index(self.jat+1) - self._T_from_index(self.jat)
            d_tin_i = 1.0 / d_tin

            psi[0] = self._psi0(y)
            psi[1] = self._psi0(1-y)
            psi[2] = self._psi1(y)
            psi[3] = self._psi1(1-y)

            din_coeff[0] = d_tin_i
            din_coeff[1] = -d_tin_i
            din_coeff[2] = 1.0
            din_coeff[3] = 1.0

            wt = psi * din_coeff
            return wt

        if sel == "biquintic":

            wt = np.zeros(6)
            coeff_derivatives = {0: _no_derivative_coeff,
                                 1: _first_derivative_coeff,
                                 2: _second_derivative_coeff}

        elif sel == 'bicubic':

            wt = np.zeros(4)
            coeff_derivatives = {0: _no_xderivative_coeff,
                                 1: _first_xderivative_coeff}
        else:

            raise NotImplementedError("There is no implementation for the selected type")

        func = coeff_derivatives.get(der)
        wt = func(self, y)

        return wt

    def _interpolate_biquintic(self, fi, wt):

        fwtr = np.zeros(6)

        fwtr[0] = fi[ 0]*wt[0] + fi[ 9]*wt[1] + fi[ 2]*wt[2] + fi[11]*wt[3] + fi[ 4]*wt[4] + fi[13]*wt[5]
        fwtr[1] = fi[18]*wt[0] + fi[27]*wt[1] + fi[20]*wt[2] + fi[29]*wt[3] + fi[22]*wt[4] + fi[31]*wt[5]
        fwtr[2] = fi[ 1]*wt[0] + fi[10]*wt[1] + fi[ 5]*wt[2] + fi[14]*wt[3] + fi[ 7]*wt[4] + fi[16]*wt[5]
        fwtr[3] = fi[19]*wt[0] + fi[28]*wt[1] + fi[23]*wt[2] + fi[32]*wt[3] + fi[25]*wt[4] + fi[34]*wt[5]
        fwtr[4] = fi[ 3]*wt[0] + fi[12]*wt[1] + fi[ 6]*wt[2] + fi[15]*wt[3] + fi[ 8]*wt[4] + fi[17]*wt[5]
        fwtr[5] = fi[21]*wt[0] + fi[30]*wt[1] + fi[24]*wt[2] + fi[33]*wt[3] + fi[26]*wt[4] + fi[35]*wt[5]

        return fwtr

    def _interpolate_bicubic(self, fi, wt):

        xfwtr = np.zeros(4)

        xfwtr[0] = fi[ 0]*wt[0] + fi[ 4]*wt[1] + fi[ 2]*wt[2] + fi[ 6]*wt[3]
        xfwtr[1] = fi[ 8]*wt[0] + fi[12]*wt[1] + fi[10]*wt[2] + fi[14]*wt[3]
        xfwtr[2] = fi[ 1]*wt[0] + fi[ 5]*wt[1] + fi[ 3]*wt[2] + fi[ 7]*wt[3]
        xfwtr[3] = fi[ 9]*wt[0] + fi[13]*wt[1] + fi[11]*wt[2] + fi[15]*wt[3]

        return xfwtr

    def _compute_free_energy(self):

        fi00 = self.f[self.iat][self.jat]
        fi01 = self.f[self.iat][self.jat+1]
        fi10 = self.f[self.iat+1][self.jat]
        fi11 = self.f[self.iat+1][self.jat+1]

        fi = np.array([fi00, fi10, fi01, fi11]).flatten()

        wd = self._read_wd(der=0)
        wt = self._read_wt(der=0)

        fwtr = self._interpolate_biquintic(fi, wt)

        free_energy = np.sum(fwtr * wd)
        return free_energy

    def _compute_pressure(self):

        fi00 = self.f[self.iat][self.jat]
        fi01 = self.f[self.iat][self.jat+1]
        fi10 = self.f[self.iat+1][self.jat]
        fi11 = self.f[self.iat+1][self.jat+1]

        fi = np.array([fi00, fi10, fi01, fi11]).flatten()

        wd = self._read_wd(der=1)
        wt = self._read_wt(der=0)

        fwtr = self._interpolate_biquintic(fi, wt)

        df_d = np.sum(fwtr * wd)

        pressure = self.din**2 * df_d
        return pressure

    def _compute_entropy(self):

        fi00 = self.f[self.iat][self.jat]
        fi01 = self.f[self.iat][self.jat+1]
        fi10 = self.f[self.iat+1][self.jat]
        fi11 = self.f[self.iat+1][self.jat+1]

        wd = self._read_wd(der=0)
        wt = self._read_wt(der=1)

        fi = np.array([fi00, fi10, fi01, fi11]).flatten()

        fwtr = self._interpolate_biquintic(fi,wt)
        df_t = np.sum(fwtr * wd)

        entropy = -df_t
        return entropy

    def _compute_energy(self):

        f = self._compute_free_energy()
        s = self._compute_entropy()

        energy = f + self.T*s
        return energy

    def _compute_dp_dr(self):

        fi00 = self.pd[self.iat][self.jat]
        fi01 = self.pd[self.iat][self.jat+1]
        fi10 = self.pd[self.iat+1][self.jat]
        fi11 = self.pd[self.iat+1][self.jat+1]

        fi = np.array([fi00, fi10, fi01, fi11]).flatten()

        wd = self._read_wd(der=0, sel='bicubic')
        wt = self._read_wt(der=0, sel='bicubic')

        fwtr = self._interpolate_bicubic(fi, wd)
        dpdr = np.sum(fwtr*wt)

        return dpdr

    def _compute_eta(self):

        fi00 = self.ef[self.iat][self.jat]
        fi01 = self.ef[self.iat][self.jat+1]
        fi10 = self.ef[self.iat+1][self.jat]
        fi11 = self.ef[self.iat+1][self.jat+1]

        fi = np.array([fi00, fi10, fi01, fi11]).flatten()

        wd = self._read_wd(der=0, sel='bicubic')
        wt = self._read_wt(der=0, sel='bicubic')

        fwtr = self._interpolate_bicubic(fi, wd)
        eta = np.sum(fwtr * wt)

        return eta

    def _compute_number_fraction(self):

        fi00 = self.xnf[self.iat][self.jat]
        fi01 = self.xnf[self.iat][self.jat+1]
        fi10 = self.xnf[self.iat+1][self.jat]
        fi11 = self.xnf[self.iat+1][self.jat+1]

        fi = np.array([fi00, fi10, fi01, fi11]).flatten()

        wd = self._read_wd(der=0, sel='bicubic')
        wt = self._read_wt(der=0, sel='bicubic')

        fwtr = self._interpolate_bicubic(fi, wd)
        xn = np.sum(fwtr * wt)

        return xn

    def update(self, state):

        self.din = state.rho * state.ye
        self.T = state.T
        self.ye = state.ye

        self.iat, self.jat = self._rho_T_indices()

        state.p += self._compute_pressure()
        state.s += self._compute_entropy()
        state.eint += self._compute_energy()
        state.dpdr += self._compute_dp_dr()
        state.eta += self._compute_eta()
        state.xn += self._compute_number_fraction()


if __name__== "__main__":

    #state = EosState(rho=1.e15, T=1.0e13, ye=1.0)
    state = EosState(rho=1.e-12, T=1.0e3, ye=1.0)
    table = HelmholtzTable()
    table.update(state)

    print(f"The pressure               : {state.p}")
    print(f"Entropy                    : {state.s}")
    print(f"Energy                     : {state.eint}")
    print(f"Pressure density derivative: {state.dpdr}")
    print(f"Chemical potential         : {state.eta}")
    print(f"Number Fraction            : {state.xn}")