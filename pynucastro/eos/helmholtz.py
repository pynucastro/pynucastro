import numpy as np

""" Here we implement the helmholtz eos. The idea is to call for two thermodynamic variables and the composition
of a reaction network."""

class HelmholtzTable:

    #The first step is to interpolate the electron_positron eos table by the use of biquintic
    #polinomials.

    def __init__(self, state, table='helm_table.dat'):

        self.din = state.rho * state.ye
        self.T = state.T

        self.imax = 541   #541 number of dens pts
        self.tlo = 3.0
        self.thi = 13.0
        self.tstp  = (self.thi - self.tlo) / (self.jmax-1)
        self.tstpi = 1.0 / self.tstp

        self.jmax = 201   #201 number of temp pts
        self.dlo   = -12.0
        self.dhi   = 15.0
        self.dstp  = (self.dhi - self.dlo) / (self.imax-1)
        self.dstpi = 1.0 / self.dstp

        self.iat = None
        self.jat = None
        self.rho_T_indices()

        self.table=table

    # Now, let us define a table row index function

    def rho_T_indices(self):

        iat = int((np.log10(self.din) - self.dlo) * self.dstpi)
        iat = max(0, min(iat, self.imax-2))

        jat = int((np.log10(self.T) - self.tlo) * self.tstpi)
        jat = max(0, min(jat, self.jmax-2))

        return iat, jat

    def indices_to_row(self, iat, jat):

        return jat + self.jmax * iat

    def row_number(self):

        iat, jat = self.rho_T_indices()
        index = self.indices_to_row(iat, jat)

        return index

    def _read_fi(self):
        """ From here we read the table only once to interpolate F"""

        def _extract_line_data(line):
            line_entries = line.strip().split()
            f = np.array([np.float64(s) for s in line_entries])
            return f

        i0, j0 = self.rho_T_indices()
        idx00 = self.indices_to_row(i0, j0)
        idx01 = self.indices_to_row(i0, j0+1)
        idx10 = self.indices_to_row(i0+1, j0)
        idx11 = self.indices_to_row(i0+1, j0+1)

        funcs: dict = {idx00: _extract_line_data,
                       idx01: _extract_line_data,
                       idx10: _extract_line_data,
                       idx11: _extract_line_data}

        with open(self.table, 'r') as file:

            f = []
            for line_counter, line in enumerate(file):

                extract_line = funcs.get(line_counter)

                if extract_line:
                    f.append(extract_line(line))

                if line_counter > idx11:
                    break

            fi = np.reshape(np.array(f), 36)

        return fi

    def _read_wt(self):

        return None

    #Let us define all the polynomials that will be defined
    #=========

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

    #=============================
    def rho_from_index(self, i):

        logrho_i = self.dlo + i * self.dstp
        return 10.0**logrho_i

    def T_from_index(self, j):

        logT_j = self.tlo + j * self.tstp
        return 10.0**logT_j

    def read_d(self, der=0):

        din_coeff = np.zeros(6)
        psi = np.zeros(6)
        x = (self.rho - self.rho_from_index(self.iat)) / \
            (self.rho_from_index(self.iat+1) - self.rho_from_index(self.iat))
        wd = np.zeros(6)

        def no_derivative_coeff(self, psi, x, din_coeff):

            wd = np.zeros(6)

            d_din = self.rho_from_index(self.iat+1) - self.rho_from_index(self.iat)
            d2_din = (self.rho_from_index(self.iat+1) - self.rho_from_index(self.iat))**2

            psi[0] = self._psi0(x)
            psi[1] = self._psi0(1-x)
            psi[2] = self._psi1(x)
            psi[3] = self._psi1(1-x)
            psi[4] = self._psi2(x)
            psi[5] = self._psi2(1-x)

            din_coeff[0] = 1.0
            din_coeff[1] = 1.0
            din_coeff[2] = d_din
            din_coeff[3] = - d_din
            din_coeff[4] = d2_din
            din_coeff[5] = d2_din

            wd = psi * din_coeff
            return wd

        def first_derivative_coeff(self, psi, x, din_coeff):

            wd = np.zeros(6)

            d_din = self.rho_from_index(self.iat+1) - self.rho_from_index(self.iat)
            d_din_i = 1.0 / d_din

            psi[0] = self._dpsi0(x)
            psi[1] = self._dpsi0(1-x)
            psi[2] = self._dpsi1(x)
            psi[3] = self._dpsi1(1-x)
            psi[4] = self._dpsi2(x)
            psi[5] = self._dpsi2(1-x)

            din_coeff[0] = d_din_i
            din_coeff[1] = - d_din_i
            din_coeff[2] = 1.0
            din_coeff[3] = 1.0
            din_coeff[4] = d_din
            din_coeff[5] = d_din

            wd = psi * din_coeff
            return wd

        def second_derivative_coeff(self, psi, x, din_coeff):

            wd = np.zeros(6)

            d_din = self.rho_from_index(self.iat+1) - self.rho_from_index(self.iat)
            d_din_i = 1.0 / d_din

            d2_din = (self.rho_from_index(self.iat+1) - self.rho_from_index(self.iat))**2
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

        coeff_derivatives = {0: no_derivative_coeff,
                             1: first_derivative_coeff,
                             2: second_derivative_coeff}

        func = coeff_derivatives.get(der)
        wd = func(psi, x, din_coeff)

        return wd

    def read_t(self, der=0):

        tin_coeff = np.zeros(6)
        psi = np.zeros(6)
        y = (self.T - self.T_from_index(self.iat)) / \
            (self.T_from_index(self.iat+1) - self.T_from_index(self.iat))
        wt = np.zeros(6)

        def no_derivative_coeff(self, psi, y, tin_coeff):

            wt = np.zeros(6)

            d_tin = self.T_from_index(self.iat+1) - self.T_from_index(self.iat)
            d2_tin = (self.T_from_index(self.iat+1) - self.T_from_index(self.iat))**2

            psi[0] = self._psi0(y)
            psi[1] = self._psi0(1-y)
            psi[2] = self._psi1(y)
            psi[3] = self._psi1(1-y)
            psi[4] = self._psi2(y)
            psi[5] = self._psi2(1-y)

            tin_coeff[0] = 1.0
            tin_coeff[1] = 1.0
            tin_coeff[2] = d_tin
            tin_coeff[3] = - d_tin
            tin_coeff[4] = d2_tin
            tin_coeff[5] = d2_tin

            wt = psi * tin_coeff
            return wt

        def first_derivative_coeff(self, psi, y, tin_coeff):

            wd = np.zeros(6)

            d_tin = self.T_from_index(self.iat+1) - self.T_from_index(self.iat)
            d_tin_i = 1.0 / d_tin

            psi[0] = self._dpsi0(y)
            psi[1] = self._dpsi0(1-y)
            psi[2] = self._dpsi1(y)
            psi[3] = self._dpsi1(1-y)
            psi[4] = self._dpsi2(y)
            psi[5] = self._dpsi2(1-y)

            tin_coeff[0] = d_tin_i
            tin_coeff[1] = - d_tin_i
            tin_coeff[2] = 1.0
            tin_coeff[3] = 1.0
            tin_coeff[4] = d_tin
            tin_coeff[5] = d_tin

            wd = psi * tin_coeff
            return wd

        def second_derivative_coeff(self, psi, y, tin_coeff):

            wd = np.zeros(6)

            d_tin = self.T_from_index(self.iat+1) - self.T_from_index(self.iat)
            d_tin_i = 1.0 / d_tin

            d2_tin = (self.T_from_index(self.iat+1) - self.T_from_index(self.iat))**2
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

            wd = psi * tin_coeff
            return wd

        coeff_derivatives = {0: no_derivative_coeff,
                             1: first_derivative_coeff,
                             2: second_derivative_coeff}

        func = coeff_derivatives.get(der)
        wt = func(psi, y, tin_coeff)

        return wt

    def interpolate_biquintic(self, fi, wt):

        fwtr = np.zeros(6)

        fwtr[0] = fi[ 0]*wt[0] + fi[ 1]*wt[1] + fi[ 2]*wt[2] + fi[ 9]*wt[3] + fi[10]*wt[4] + fi[11]*wt[5]
        fwtr[1] = fi[18]*wt[0] + fi[19]*wt[1] + fi[20]*wt[2] + fi[27]*wt[3] + fi[28]*wt[4] + fi[29]*wt[5]
        fwtr[2] = fi[ 3]*wt[0] + fi[ 5]*wt[1] + fi[ 7]*wt[2] + fi[12]*wt[3] + fi[14]*wt[4] + fi[16]*wt[5]
        fwtr[3] = fi[21]*wt[0] + fi[23]*wt[1] + fi[25]*wt[2] + fi[30]*wt[3] + fi[32]*wt[4] + fi[34]*wt[5]
        fwtr[4] = fi[ 4]*wt[0] + fi[ 6]*wt[1] + fi[ 8]*wt[2] + fi[13]*wt[3] + fi[15]*wt[4] + fi[17]*wt[5]
        fwtr[5] = fi[22]*wt[0] + fi[24]*wt[1] + fi[26]*wt[2] + fi[31]*wt[3] + fi[33]*wt[4] + fi[35]*wt[5]

        return fwtr

    def compute_free_energy(self):

        fi = self._read_fi()
        wd = self.read_d(der=0)
        wt = self.read_t(der=0)

        fwtr = self.interpolate_biquintic(fi,wt)

        free_energy = np.sum(fwtr * wd)
        return free_energy

    def compute_pressure(self):

        fi = self._read_fi()
        wd = self.read_d(der=1)
        wt = self.read_t(der=0)

        fwtr = self.interpolate_biquintic(fi,wt)
        dp_d = np.sum(fwtr * wd)

        pressure = self.din**2 * dp_d

        return pressure

    def compute_entropy(self):

        fi = self._read_fi()
        wd = self.read_d(der=0)
        wt = self.read_t(der=1)

        fwtr = self.interpolate_biquintic(fi,wt)
        df_t = np.sum(fwtr * wd)

        entropy = -df_t

        return entropy

    def compute_energy(self):

        fi = self._read_fi()
        wd = self.read_d(der=1)
        wt = self.read_t(der=1)

        fwtr = self.interpolate_biquintic(fi,wt)

        f = self.compute_free_energy()
        s = self.compute_entropy

        energy = f + self.T*s
        return energy

    def compute_chemical_potential(self):

        raise NotImplementedError








    def interpolate_bitric(self, fi, wt):

        raise NotImplementedError


if __name__== "__main__":

    # table = HelmholtzTable(rho=1.e15, T=1.0e13, ye=1.0)
    # n = table.row_number()
    # i, j = table.rho_T_indices()
    # print(i,j)
    # print(n)

    # table = HelmholtzTable(rho=1.e-12, T=1.0e3, ye=1.0)
    # n = table.row_number()
    # i, j = table.rho_T_indices()
    # print(i,j)
    # print(n)

    #table = HelmholtzTable(rho=1.e15, T=1.0e13, ye=1.0)
    # table = HelmholtzTable(rho=1.e-12, T=1.0e3, ye=1.0)
    # f = table._read_f()
    # print(f)

    pass