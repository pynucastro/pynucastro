"""Classes and methods that provide access to the AME mass excess data."""

from pathlib import Path


class MassTable:
    """Store the nuclear masses as a dictionary mapping nuclides to
    their mass excess (A_nuc - A measured) in MeV from the table
    mass_excess2020.txt.

    Parameters
    ----------
    filename : str, pathlib.Path
        The name of the file containing nuclei and mass excesses
        (mass_excess2020.txt is used by default).

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
        """Return the mass excess for a nucleus

        Parameters
        ----------
        a : int
            Atomic weight
        z : int
            Atomic number

        Returns
        -------
        float

        """

        try:
            return self.mass_diff[a, z]
        except KeyError as exc:
            raise NotImplementedError(f"nuclear mass difference for A={a} and Z={z} not available") from exc
