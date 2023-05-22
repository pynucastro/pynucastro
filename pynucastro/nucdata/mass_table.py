import os


class MassTable:
    """
    The purpose of this class is to:

    a) Create a dictionary mapping nuclides to their mass excess A_nuc - A measured
       in MeV from the table mass_excess2020.txt.

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
            datafile_name = 'mass_excess2020.txt'
            nucdata_dir = os.path.dirname(os.path.realpath(__file__))
            self.filename = os.path.join(os.path.join(nucdata_dir, 'AtomicMassEvaluation'), datafile_name)

        self._read_table()

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

            self._mass_diff[A, Z] = dm

        file.close()

    def get_mass_diff(self, a, z):
        if (a, z) in self._mass_diff:
            return self._mass_diff[a, z]
        raise NotImplementedError("Nuclear mass difference is not available")
