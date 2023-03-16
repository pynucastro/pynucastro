import os


class SpinTable:
    """
    This class stores the nubase2020_1.txt table information in a dictionary data structure
    that maps nuclides to their number of nuclear spin states.
    Therefore, after setting an SpinTable class in rates.py, we can retrieve the
    spin states for a designated Nucleus class.

    The variable reliable switch between using all the values of the tables, excluding the nuclei
    where only intervals are given and the values measured by strong experimental arguments.
    """

    def __init__(self, datafile=None, reliable=False):

        self._spin_states = {}
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
                    self._spin_states[A, Z] = spin_states
                else:
                    continue
            else:
                self._spin_states[A, Z] = spin_states

        finput.close()

    def get_spin_states(self, a, z):
        if (a, z) in self._spin_states:
            return self._spin_states[a, z]
        raise NotImplementedError("nuclear spin data is not available for the selected nucleus")
