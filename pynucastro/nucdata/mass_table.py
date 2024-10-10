from pathlib import Path


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

    def __init__(self, filename: str | Path = None):

        self.mass_diff = {}

        if filename:
            self.filename = Path(filename)
        else:
            nucdata_dir = Path(__file__).parent
            datafile_name = 'mass_excess2020.txt'
            self.filename = nucdata_dir/'AtomicMassEvaluation'/datafile_name

        self._read_table()

    def _read_table(self) -> None:

        with self.filename.open("r") as f:
            # skip the header
            for line in f.readlines()[5:]:
                A, Z, dm = line.strip().split()[:3]
                self.mass_diff[int(A), int(Z)] = float(dm)

    def get_mass_diff(self, a: int, z: int) -> float:
        try:
            return self.mass_diff[a, z]
        except KeyError as exc:
            raise NotImplementedError(f"nuclear mass difference for A={a} and Z={z} not available") from exc
