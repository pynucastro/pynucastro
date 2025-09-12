"""Classes and methods that provide access to the AME half life data."""

from pathlib import Path


class HalfLifeTable:
    """Read the table of halflives (in seconds) extracted from the
    Nubase compilation in halflife2020.txt

    """

    def __init__(self, filename: str | Path = None) -> None:

        self.halflife = {}

        if filename:
            self.filename = Path(filename)
        else:
            nucdata_dir = Path(__file__).parent
            datafile_name = 'halflife2020.txt'
            self.filename = nucdata_dir/'AtomicMassEvaluation'/datafile_name

        self._read_table()

    def _read_table(self) -> None:

        with open(self.filename, "r") as f:
            # skip the header
            for line in f.readlines()[5:]:

                A, Z, tau = line.strip().split()[:3]
                #print(data_list)

                if tau != "stable":
                    tau = float(tau)

                self.halflife[int(A), int(Z)] = tau

    def get_halflife(self, a: int, z: int) -> float | str:
        """Return the half life of a nucleus.

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
            return self.halflife[a, z]
        except KeyError as exc:
            raise NotImplementedError(f"Nuclear halflife for A={a} and Z={z} not available") from exc
