from pynucastro.nucdata import PeriodicTable
import os

filename = 'mass_excess2020.txt'
dir_nucdata = os.path.dirname(os.path.realpath(__file__))
dir_mass_data = os.path.join(os.path.join(dir_nucdata, 'AtomicMassEvaluation'), filename)


class MassNuclide:

    def __init__(self, a, z, dm):
        self.a = a
        self.z = z
        self.dm = dm

    def __repr__(self):

        if self.a == 1 and self.z == 0:
            rep = 'n'
        else:
            element = PeriodicTable.lookup_Z(self.z)
            rep = '{}{}'.format(element.abbreviation, self.a)

        return rep

    def __eq__(self, other):

        return self.a == other.a and self.z == self.z


class MassTable:

    def __init__(self, filename=None):

        self._mass_diff = {}

        if filename:
            self.filename = filename
        else:
            self.filename = dir_mass_data

        self._read_table()

    def _add_mass_nuclide(self, nuc):

        assert isinstance(nuc, MassNuclide)
        assert not str(nuc) in self._mass_diff.keys()

        self._mass_diff[str(nuc)] = nuc

    def _read_table(self):

        file = open(self.filename, 'r')

        for _ in range(4):
            file.readline()

        for line in file:

            data_list = line.strip().split()
            #print(data_list)
            A_str = data_list.pop(0)
            Z_str = data_list.pop(0)
            dm_str = data_list.pop(0)

            A = int(A_str)
            Z = int(Z_str)
            dm = float(dm_str)

            nuc = MassNuclide(a=A, z=Z, dm=dm)
            self._add_mass_nuclide(nuc)

        file.close()

    def get_mass_diff(self, nuc):

        if str(nuc) in self._mass_diff.keys():
            return self._mass_diff[str(nuc)]
        else:
            raise NotImplementedError("Nuclear mass difference is not available")
