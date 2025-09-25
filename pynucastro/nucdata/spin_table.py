"""Classes and methods for dealing with nuclear spin data."""

from pathlib import Path


class SpinTable:
    """Manage the nuclear spin data in a dictionary data structure
    that maps nuclides to their number of nuclear spin states.
    Therefore, after setting an SpinTable class in rates.py, we can
    retrieve the spin states for a designated Nucleus class.

    The variable reliable switch between using all the values of the
    tables, excluding the nuclei where only intervals are given and
    the values measured by strong experimental arguments.

    """

    def __init__(self, datafile: str | Path = None) -> None:

        self._spin_states = {}
        self._reliability_table = {}

        if datafile:
            self.datafile = Path(datafile)
        else:
            nucdata_dir = Path(__file__).parent
            datafile_name = 'spins2020.txt'
            self.datafile = nucdata_dir/'AtomicMassEvaluation'/datafile_name

        self._read_table()

    def _read_table(self) -> None:

        with self.datafile.open("r") as f:

            for line in f.readlines()[4:]:

                A, Z, _, spin_states, experimental = line.strip().split()[:5]
                A, Z, spin_states = int(A), int(Z), int(spin_states)

                self._spin_states[A, Z] = spin_states
                self._reliability_table[A, Z] = experimental == 's'

    def get_spin_states(self, a: int, z: int, reliable: bool = False) -> int:
        """Return the spin for a nucleus.

        Parameters
        ----------
        a : int
            Atomic weight
        z : int
            Atomic number
        reliable : bool
            setting this to True will only return spin states
            for nuclei where the spin is known from strong
            experimental arguments

        Returns
        -------
        float

        """
        if reliable:
            try:
                if self._reliability_table[a, z]:
                    return self._spin_states[a, z]
                else:
                    raise NotImplementedError(f"nuclear spin data for A={a} and Z={z} not reliable enough for your specifications")
            except KeyError as exc:
                raise NotImplementedError(f"nuclear spin data for A={a} and Z={z} not available") from exc

        try:
            return self._spin_states[a, z]
        except KeyError as exc:
            raise NotImplementedError(f"nuclear spin data for A={a} and Z={z} not available") from exc
