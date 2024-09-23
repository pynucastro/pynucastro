import os


class HalfLifeTable:
    """
    Read the table of halflives (in seconds) extracted from the Nubase
    compilation in halflife2020.txt
    """

    def __init__(self, filename=None):

        self.halflife = {}

        if filename:
            self.filename = filename
        else:
            datafile_name = 'halflife2020.txt'
            nucdata_dir = os.path.dirname(os.path.realpath(__file__))
            self.filename = os.path.join(os.path.join(nucdata_dir, 'AtomicMassEvaluation'), datafile_name)

        self._read_table()

    def _read_table(self):

        file = open(self.filename, 'r')

        # skip the header
        for _ in range(5):
            file.readline()

        for line in file:

            data_list = line.strip().split()
            #print(data_list)
            A_str = data_list.pop(0)
            Z_str = data_list.pop(0)
            tau = data_list.pop(0)

            A = int(A_str)
            Z = int(Z_str)
            if tau != "stable":
                tau = float(tau)

            self.halflife[A, Z] = tau

        file.close()

    def get_halflife(self, a, z):
        if (a, z) in self.halflife:
            return self.halflife[a, z]
        raise NotImplementedError("Nuclear halflife is not available")
