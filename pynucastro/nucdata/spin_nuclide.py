import os
import numpy as np

from pynucastro.nucdata import PeriodicTable


class SpinNuclide:
    """
    This class contains the spin information for each nuclei specified in the nubase20.txt table.
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
    This class stores the nubase2020.txt table information in a dictionary data structure
    that maps a pynucastro.rates.Nucleus object to a pynucastro.nucdata.SpinNuclide object.
    Therefore, after setting an SpinTable class in rates.py, we should retrieve the SpinNuclide
    data structure from a designated Nucleus class.
    """

    def __init__(self, datafile=None, set_double_gs=False):

        self.set_double_gs = set_double_gs
        self._spin_nuclide = {}

        if datafile:
            self.datafile = datafile
        else:
            datafile_name = 'nubase2020.txt'
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
            state_list = []

            if len(ls) == 6:
                if not self.set_double_gs:
                    continue
                A = int(ls.pop(0))
                Z = int(ls.pop(0))
                ls.pop(0)           # Spin 1
                ls.pop(0)           # Spin 2
                states1 = ls.pop(0)
                states2 = ls.pop(0)

                state_list.append(states1)
                state_list.append(states2)
                double_spin_states = np.array([float(s) for s in state_list])
                spin_nuc = SpinNuclide(a=A, z=Z, spin_states=double_spin_states)
                self._add_spin_nuclide(spin_nuc)
            else:
                A = int(ls.pop(0))
                Z = int(ls.pop(0))
                ls.pop(0)           # Spin
                spin_states = float(ls.pop(0))
                spin_nuc = SpinNuclide(a=A, z=Z, spin_states=spin_states)
                self._add_spin_nuclide(spin_nuc)

        finput.close()

    def get_spin_nuclide(self, nuc):
        if str(nuc) in self._spin_nuclide:
            return self._spin_nuclide[str(nuc)]
        raise NotImplementedError("nuclear spin data is not available for the selected nucleus")
