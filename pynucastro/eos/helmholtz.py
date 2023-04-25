import numpy as np

""" Here we implement the helmholtz eos. The idea is to call for two thermodynamic variables and the composition
of a reaction network."""

class HelmholtzTable:

    #The first step is to interpolate the electron_positron eos table by the use of biquintic
    #polinomials.

    def __init__(self, rho, T, ye, table='helm_table.dat'):

        self.rho = rho
        self.T = T
        self.ye = ye

        self.imax = 541   #541 number of dens pts
        self.jmax = 201   #201 number of temp pts
        self.tlo = 3.0
        self.thi = 13.0
        self.tstp  = (self.thi - self.tlo) / (self.jmax-1)
        self.tstpi = 1.0 / self.tstp
        self.dlo   = -12.0
        self.dhi   = 15.0
        dstp  = (self.dhi - self.dlo) / (self.imax-1)
        self.dstpi = 1.0 / dstp

        self.iat = None
        self.jat = None

        self.table=table

    # Now, let us define a table row index function

    def rho_T_indices(self):

        din = self.rho * self.ye

        iat = int((np.log10(din) - self.dlo) * self.dstpi)
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

    def interpolate_biquintic(self, fi, wt):

        #fi = self._read_fi()
        #sit, sid = self._read_td()

        fwtr = np.zeros(6)

        fwtr[0] = fi[] * wt[0] + fi[] * wt[1] + fi[] * wt[2] + fi[] * wt[3] + fi * wt[4] + fi[] * wt[5]
        fwtr[1] = fi[] * wt[0] + fi[] * wt[1] + fi[] * wt[2] + fi[] * wt[3] + fi * wt[4] + fi[] * wt[5]
        fwtr[2] = fi[] * wt[0] + fi[] * wt[1] + fi[] * wt[2] + fi[] * wt[3] + fi * wt[4] + fi[] * wt[5]
        fwtr[3] = fi[] * wt[0] + fi[] * wt[1] + fi[] * wt[2] + fi[] * wt[3] + fi * wt[4] + fi[] * wt[5]
        fwtr[4] = fi[] * wt[0] + fi[] * wt[1] + fi[] * wt[2] + fi[] * wt[3] + fi * wt[4] + fi[] * wt[5]
        fwtr[5] = fi[] * wt[0] + fi[] * wt[1] + fi[] * wt[2] + fi[] * wt[3] + fi * wt[4] + fi[] * wt[5]

        return fwtr

    def interpolate_bitric(self, fi, wt):

        raise NotImplementedError
