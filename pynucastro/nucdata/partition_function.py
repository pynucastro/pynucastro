"""Classes and methods for dealing with nuclear partition functions."""


from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline


class PartitionFunction:
    """Store the tabulated data for the partition function for a
    specific nucleus, which can be combined with other
    (non-overlapping) partition functions by addition and evaluated
    for different temperature values.

    Adding two PartitionFunction objects is implemented by simply appending the
    temperature and partition function arrays of the higher-temperature
    partition function to those of the lower-temperature partition function. If
    the temperature ranges overlap, however, an exception is generated.

    If either of the PartitionFunction objects added have already had a spline
    interpolant constructed, then construct a new spline interpolant for the
    returned PartitionFunction of order equal to the maximum order of the added
    PartitionFunction objects.

    Parameters
    ----------
    nucleus : str
        The nucleus (e.g. ``"ni56"``)
    name : str
        The name of the table on which the nucleus is read
    temperature : numpy.ndarray
        A sorted array of all the temperatures involved
    partition_function : numpy.ndarray
        An array with all the partition function values given in the
        same order as ``temperature``
    interpolant_order : int
        The interpolation spline order, must be between 1 and 5, inclusive

    """

    def __init__(self, nucleus, name, temperature,
                 partition_function, interpolant_order=3):
        assert isinstance(nucleus, str)

        temperature = np.asarray(temperature)
        partition_function = np.asarray(partition_function)

        assert temperature.shape == partition_function.shape
        assert np.all(temperature[:-1] <= temperature[1:]), "temperature array must be sorted"

        self.nucleus = nucleus
        self.name = name
        self.temperature = temperature
        self.partition_function = partition_function
        self.interpolant_order = interpolant_order
        self._interpolant = None

    def lower_temperature(self):
        """Return the lowest temperature this object supports.

        Returns
        -------
        float

        """

        return self.temperature[0]

    def upper_temperature(self):
        """Return the highest temperature this object supports.

        Returns
        -------
        float
        """

        return self.temperature[-1]

    def __add__(self, other):
        assert self.nucleus == other.nucleus

        if self.upper_temperature() < other.lower_temperature():
            lower = self
            upper = other
        else:
            lower = other
            upper = self
        if lower.upper_temperature() >= upper.lower_temperature():
            raise ValueError("temperature ranges cannot overlap")

        temperature = np.concatenate([lower.temperature, upper.temperature])
        partition_function = np.concatenate([lower.partition_function,
                                             upper.partition_function])

        name = f'{lower.name}+{upper.name}'

        order = max(self.interpolant_order, other.interpolant_order)
        newpf = PartitionFunction(nucleus=self.nucleus, name=name,
                                  temperature=temperature,
                                  partition_function=partition_function,
                                  interpolant_order=order)

        return newpf

    def __eq__(self, other):
        return (np.all(self.partition_function == other.partition_function) and
                np.all(self.temperature == other.temperature))

    def eval(self, T):
        """Compute the interpolated partition function value for the
        temperature T.

        Parameters
        ----------
        T : float
            temperature in K

        Returns
        -------
        float

        """

        # lazily construct the interpolant object, since it's pretty expensive
        if not self._interpolant:
            self._interpolant = InterpolatedUnivariateSpline(
                self.temperature/1.0e9,
                np.log10(self.partition_function),
                k=self.interpolant_order
            )
        try:
            T = float(T)/1.0e9
        except ValueError:
            print("invalid temperature")
            raise
        # extrapolates keeping the boundaries fixed.
        return float(10**self._interpolant(T, ext='const'))

    def plot(self, T_min=1.e7, T_max=1.e11):
        """Plot the partition function as a function of T

        Parameters
        ----------
        T_min : float
            minimum temperature to plot
        T_max : float
            maximum temperature to plot

        Returns
        -------
        matplotlib.figure.Figure

        """

        fig, ax = plt.subplots()

        ax.plot(self.temperature, self.partition_function)
        mask = (self.temperature >= T_min) & (self.temperature <= T_max)
        ax.grid(ls=":")
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlim(T_min, T_max)
        ax.set_ylim(min(0.5, self.partition_function[mask].min()),
                    max(2, self.partition_function[mask].max()))
        ax.set_xlabel("T (K)")
        ax.set_ylabel("G(T)")

        return fig

class PartitionFunctionTable:
    """Manage a partition function table file. A
    :class:`PartitionFunction` object is constructed for each nucleus
    and stored in a dictionary keyed by the lowercase nucleus name in
    the form e.g.  "ni56". The table files are stored in the
    ``PartitionFunction`` subdirectory.

    :var name:         the name of the table (as defined in the data file)
    :var temperatures: an array of temperature values

    """

    def __init__(self, file_name):
        self._partition_function = {}
        self.name = None
        self.temperatures = None
        self._read_table(file_name)

    def _add_nuclide_pfun(self, nuc, pfun):
        assert isinstance(nuc, str)
        assert nuc not in self._partition_function
        self._partition_function[nuc] = pfun

    def get_nuclei(self):
        """Return a set of the nuclei this table supports."""
        return set(self._partition_function)

    def get_partition_function(self, nuc):
        """Return the :class:`PartitionFunction` object for a specific nucleus."""
        assert isinstance(nuc, str)
        if nuc in self._partition_function:
            return self._partition_function[nuc]
        return None

    def _read_table(self, file_name: str | Path):
        with Path(file_name).open("r") as fin:

            # get headers name
            fhead = fin.readline()
            hsplit = fhead.split('name: ')
            self.name = hsplit[-1].strip('\n')

            # throw away the six subsequent lines
            for _ in range(6):
                fin.readline()

            # Now, we want to read the lines of the file where
            # the temperatures are located
            temp_strings = fin.readline().strip().split()
            self.temperatures = np.array(temp_strings, dtype=np.float64)

            # Now, we append on the array lines = [] all the remaining file, the structure
            # 1. The nucleus
            # 2. The partition value relative to the nucleus defined in 1.

            lines = []
            for line in fin:
                ls = line.strip()
                if ls:
                    lines.append(ls)

        # Using .pop(0) twice we construct each nucleus partition function.
        while lines:
            nuc = lines.pop(0)
            pfun_strings = lines.pop(0).split()
            partitionfun = np.array(pfun_strings, dtype=np.float64)
            pfun = PartitionFunction(nuc, self.name, self.temperatures, partitionfun)
            self._add_nuclide_pfun(nuc, pfun)


class PartitionFunctionCollection:
    """A collection of :class:`PartitionFunctionTable` objects in a
    dictionary keyed by the name of the tables.

    In our discussion we have two different sets of tables: FRDM and ETFSI-Q.

    :var use_high_temperatures: whether to incorporate the high-temperature data
                                tables
    :var use_set: selects between the FRDM (``'frdm'``) and ETFSI-Q
                  (``'etfsiq'``) data sets.

    """

    def __init__(self, use_high_temperatures=True, use_set='frdm'):
        self._partition_function_tables = {}
        self.use_high_temperatures = use_high_temperatures
        self.use_set = use_set
        self._read_collection()

    def _add_table(self, table):
        """Add a PartitionFunctionTable to this collection."""
        assert table.name not in self._partition_function_tables
        self._partition_function_tables[table.name] = table

    def _read_collection(self):
        """Read and construct all the tables from the data files."""

        nucdata_dir = Path(__file__).parent
        partition_function_dir = nucdata_dir/'PartitionFunction'

        pft = PartitionFunctionTable(file_name=partition_function_dir/'etfsiq_low.txt')
        self._add_table(pft)

        pft = PartitionFunctionTable(file_name=partition_function_dir/'frdm_low.txt')
        self._add_table(pft)

        pft = PartitionFunctionTable(file_name=partition_function_dir/'etfsiq_high.txt')
        self._add_table(pft)

        pft = PartitionFunctionTable(file_name=partition_function_dir/'frdm_high.txt')
        self._add_table(pft)

    def get_nuclei(self):
        """Return a set of all the nuclei this collection supports."""
        nuclei = set()
        for table in self._partition_function_tables.values():
            nuclei.update(table.get_nuclei())
        return nuclei

    def __iter__(self):
        for nuc in self.get_nuclei():
            yield self.get_partition_function(nuc)

    def get_partition_function(self, nuc):
        """Return the :class:`PartitionFunction` object for a specific nucleus."""
        assert isinstance(nuc, str)

        if self.use_set == 'frdm':
            pf_lo_table = self._partition_function_tables['frdm_low']
            pf_lo = pf_lo_table.get_partition_function(nuc)

            pf_hi_table = self._partition_function_tables['frdm_high']
            pf_hi = pf_hi_table.get_partition_function(nuc)

        elif self.use_set == 'etfsiq':
            pf_lo_table = self._partition_function_tables['etfsiq_low']
            pf_lo = pf_lo_table.get_partition_function(nuc)

            pf_hi_table = self._partition_function_tables['etfsiq_high']
            pf_hi = pf_hi_table.get_partition_function(nuc)
        else:
            raise ValueError("invalid partition function type")

        if self.use_high_temperatures:
            if pf_lo and pf_hi:
                pf = pf_lo + pf_hi
            elif pf_lo:
                pf = pf_lo
            elif pf_hi:
                pf = pf_hi
            else:
                #name = 'default'
                #pf_default = PartitionFunction(nuc, name, pf_lo_table.temperatures, np.ones_like(pf_lo_table.temperatures))
                #pf = pf_default
                raise ValueError
                #if str(nuc) != 'h1' and str(nuc) != 'n' and str(nuc) != 'he4':
                #    print(f'WARNING: {nuc} partition function is not supported: set pf = 1.0 by default')

        else:
            if pf_lo:
                pf = pf_lo
            else:
                raise ValueError
                #name = 'default'
                #pf_default = PartitionFunction(nuc, name, pf_lo_table.temperatures, np.ones_like(pf_lo_table.temperatures))
                #pf = pf_default
                #if str(nuc) != 'p' and str(nuc) != 'n' and str(nuc) != 'he4':
                #    print(f'WARNING: {nuc} partition function is not supported: set pf = 1.0 by default')

        return pf
