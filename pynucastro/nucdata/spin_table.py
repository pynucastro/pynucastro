import os


class SpinData:
    """
    This class contains the relevant spin information for a nuclide specified in
    the nubase2020_1.txt table.

    Currently, we only use the number of spin states.
    """

    def __init__(self, spin_states):
        self.spin_states = spin_states


class SpinTable:
    """
    This class stores the nubase2020_1.txt table information in a dictionary data structure
    that maps nuclides to a pynucastro.nucdata.SpinData object.
    Therefore, after setting an SpinTable class in rates.py, we should retrieve the SpinData
    data structure for a designated Nucleus class.

    The variable reliable switch between using all the values of the tables, excluding the nuclei
    where only intervals are given and the values measured by strong experimental arguments.
    """

    def __init__(self, datafile=None, reliable=False):

        self._spin_data = {}
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

    def _read_table(self):

        finput = open(self.datafile, 'r')

        for _ in range(4):
            finput.readline()

        for line in finput:

            ls = line.strip().split()

            A = int(ls.pop(0))
            Z = int(ls.pop(0))
            ls.pop(0)           # Spin
            spin_states = int(ls.pop(0))
            experimental = ls.pop(0)

            if self.reliable:
                if experimental == 's':
                    self._spin_data[A, Z] = SpinData(spin_states)
                else:
                    continue
            else:
                self._spin_data[A, Z] = SpinData(spin_states)

        finput.close()

    def get_spin_data(self, a, z):
        if (a, z) in self._spin_data:
            return self._spin_data[a, z]
        raise NotImplementedError("nuclear spin data is not available for the selected nucleus")
