from pynucastro.nucdata import PeriodicTable
import os

filename = 'mass_excess2020.txt'
dir_nucdata = os.path.dirname(os.path.realpath(__file__))
dir_mass_data = os.path.join(os.path.join(dir_nucdata, 'AtomicMassEvaluation'), filename)


class MassNuclide:
    """
    The purpose of this class is to contain the nuclei information and create an instance
    represented by a string that contains the nuclide and the atomic weight A number, for
    example, o16, to the Oxygen-16 nuclei.

    The variables required to construct MassNuclide class instance are:

    :var a:        The Atomic weight, measured in atomic number units.
    :var z:        The Atomic Number.
    :var dm:       The mass difference, measured in MeV.
    """

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
    """
    The purpose of this class is to:

    a) Create a dictionary between MassNuclide objects, previously constructed
       from the table mass_excess2020.txt, and their mass excess A_nuc - A
       measured in MeV.

    b) Parse the information of the previously defined dictionary to the Nucleus
       class.

    The only required variable to construct an instance of this class is : var filename:
    that contains the .txt table file with the nuclei and their mass excess. If this
    variable is not provided, then mass_excess2020.txt is considered by default.
    """

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
