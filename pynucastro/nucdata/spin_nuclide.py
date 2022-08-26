import os

from pynucastro.nucdata.elements import PeriodicTable


class SpinNuclide:
    """
    This class contains the spin information for each nuclei specified in the nubase2020_1.txt table.
    The purpose to create this class lies on the
    """

    def __init__(self, a=None, z=None, spin_states=None):

        self.a = int(a)
        self.z = int(z)
        self.spin_states = spin_states

    def __repr__(self):

        if self.a == 1 and self.z == 0:
            rep = 'n'
        else:
            element = PeriodicTable.lookup_Z(self.z)
            rep = f'{element.abbreviation}{self.a}'

        return rep

    def __eq__(self, other):

        return self.a == other.a and self.z == other.z


class SpinTable:
    """
    This class stores the nubase2020_1.txt table information in a dictionary data structure
    that maps a pynucastro.rates.Nucleus object to a pynucastro.nucdata.SpinNuclide object.
    Therefore, after setting an SpinTable class in rates.py, we should retrieve the SpinNuclide
    data structure from a designated Nucleus class.

    The variable reliable switch between using all the values of the tables, excluding the nuclei
    where only intervals are given and the values measured by strong experimental arguments.
    """

    def __init__(self, datafile=None, reliable=None):

        self._spin_nuclide = {}
        self.reliable = reliable

        if datafile:
            self.datafile = datafile
        else:
            datafile_name = 'nubase2020_1.txt'
            nucdata_dir = os.path.dirname(os.path.realpath(__file__))
            datafile_dir = os.path.join(os.path.join(nucdata_dir, 'AtomicMassEvaluation'), datafile_name)

        if os.path.isfile(datafile_dir):
            self.datafile = datafile_dir
        else:
            raise Exception('ERROR: The spin tabulated file was not found')

        self._read_table()

    def _add_spin_nuclide(self, spin_nuc):

        assert isinstance(spin_nuc, SpinNuclide)
        assert str(spin_nuc) not in self._spin_nuclide

        self._spin_nuclide[str(spin_nuc)] = spin_nuc

    def _read_table(self):

        finput = open(self.datafile, 'r')

        for _ in range(4):
            finput.readline()

        for line in finput:

            ls = line.strip().split()

            A = int(ls.pop(0))
            Z = int(ls.pop(0))
            ls.pop(0)           # Spin
            spin_states = float(ls.pop(0))
            experimental = ls.pop(0)

            if self.reliable:
                if experimental == 's':
                    spin_nuc = SpinNuclide(a=A, z=Z, spin_states=spin_states)
                    self._add_spin_nuclide(spin_nuc)
                else:
                    continue
            else:
                spin_nuc = SpinNuclide(a=A, z=Z, spin_states=spin_states)
                self._add_spin_nuclide(spin_nuc)

        finput.close()

    def get_spin_nuclide(self, nuc):
        if str(nuc) in self._spin_nuclide:
            return self._spin_nuclide[str(nuc)]
        raise NotImplementedError("nuclear spin data is not available for the selected nucleus")
