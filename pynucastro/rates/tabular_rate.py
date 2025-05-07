import io
import math
from enum import Enum
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

import pynucastro.numba_util as numba
from pynucastro.nucdata import Nucleus, UnsupportedNucleus
from pynucastro.numba_util import jitclass
from pynucastro.rates.files import RateFileError, _find_rate_file
from pynucastro.rates.rate import Rate, RateSource


class TableIndex(Enum):
    """An enum-like container for indexing the electron-capture
    tables

    """
    RHOY = 0
    T = 1
    MU = 2
    DQ = 3
    VS = 4
    RATE = 5
    NU = 6
    GAMMA = 7


@jitclass([
    ('data', numba.float64[:, :]),
    ('table_rhoy_lines', numba.int32),
    ('table_temp_lines', numba.int32),
    ('rhoy', numba.float64[:]),
    ('temp', numba.float64[:])
])
class TableInterpolator:
    """A class that holds a pointer to the table data and
    methods that allow us to interpolate a variable

    Parameters
    ----------
    table_rhoy_lines : numpy.ndarray
        a 1D array giving the (ρ Y_e) values where the rate is tabulated
    table_temp_lines : numpy.ndarray
        a 1D array guving the T values where the rate is tabulated
    table_data : numpy.ndarray
        a 2D array giving the tabulated rate data of the form (index,
        component) where index is a 1D flattened representation of
        (rhoY, T).

    """

    def __init__(self, table_rhoy_lines, table_temp_lines, table_data):

        self.data = table_data
        self.table_rhoy_lines = table_rhoy_lines
        self.table_temp_lines = table_temp_lines

        # for easy indexing, store a 1-d array of T and rhoy
        self.rhoy = self.data[::self.table_temp_lines, TableIndex.RHOY.value]
        self.temp = self.data[0:self.table_temp_lines, TableIndex.T.value]

    def _get_logT_idx(self, logt0):
        """Find the index into the temperatures such that T[i-1] < t0 <= T[i].
        We return i-1 here, corresponding to the lower value.

        Note: we work in terms of log10()

        Parameters
        ----------
        logt0 : float
            log10(temperature) to interpolate at

        Returns
        -------
        int

        """

        max_idx = len(self.temp) - 1
        return max(0, min(max_idx, np.searchsorted(self.temp, logt0)) - 1)

    def _get_logrhoy_idx(self, logrhoy0):
        """Feturn the index into rhoY such that rhoY[i-1] < rhoy0 <= rhoY[i].
        We return i-1 here, corresponding to the lower value.

        Note: we work in terms of log10()

        Parameters
        ----------
        logrhoy0 : float
            log10(ρ Y_e) to interpolate at

        Returns
        -------
        int

        """

        max_idx = len(self.rhoy) - 1
        return max(0, min(max_idx, np.searchsorted(self.rhoy, logrhoy0)) - 1)

    def _rhoy_T_to_idx(self, irhoy, jtemp):
        """Given a pair (irhoy, jtemp) into the table, return the 1D
        index into the underlying data array assuming row-major
        ordering

        Parameters
        ----------
        irhoy : int
            the index in the (ρ Y_e) dimension
        ijtemp : int
            the index into the T dimension

        Returns
        -------
        int

        """

        return irhoy * self.table_temp_lines + jtemp

    def interpolate(self, logrhoy, logT, component):
        """Given logrhoy and logT, do bilinear interpolation to
        find the value of the data component in the table

        Parameters
        ----------
        logrhoy : float
            log10(ρ Y_e) to interpolate at
        logT : float
            log10(T) to interpolate at
        component : int
            the component from the data table we are interpolating.
            This should correspond to a :py:class:`TableIndex` component.

        Returns
        -------
        float

        """

        # We are going to do bilinear interpolation.  We create a
        # polynomial of the form:
        #
        # f = A [log(rho) - log(rho_i)] [log(T) - log(T_j)] +
        #     B [log(rho) - log(rho_i)] +
        #     C [log(T) - log(T_j)] +
        #     D
        #
        # we then find the i,j such that our point is in the
        # box with corners (i,j) to (i+1,j+1), and solve for
        # A, B, C, D

        # find the T and rhoY in the data table corresponding to the
        # lower left

        if logT < self.temp.min() or logT > self.temp.max():
            raise ValueError("temperature out of table bounds")

        if logrhoy < self.rhoy.min() or logrhoy > self.rhoy.max():
            raise ValueError("rhoy out of table bounds")

        irhoy = self._get_logrhoy_idx(logrhoy)
        jT = self._get_logT_idx(logT)

        # note: rhoy and T are already stored as log

        dlogrho = self.rhoy[irhoy+1] - self.rhoy[irhoy]
        dlogT = self.temp[jT+1] - self.temp[jT]

        # get the data at the 4 points

        idx = self._rhoy_T_to_idx(irhoy, jT)
        f_ij = self.data[idx, component]

        idx = self._rhoy_T_to_idx(irhoy+1, jT)
        f_ip1j = self.data[idx, component]

        idx = self._rhoy_T_to_idx(irhoy, jT+1)
        f_ijp1 = self.data[idx, component]

        idx = self._rhoy_T_to_idx(irhoy+1, jT+1)
        f_ip1jp1 = self.data[idx, component]

        D = f_ij
        C = (f_ijp1 - f_ij) / dlogT
        B = (f_ip1j - f_ij) / dlogrho
        A = (f_ip1jp1 - B * dlogrho - C * dlogT - D) / (dlogrho * dlogT)

        r = (A * (logrhoy - self.rhoy[irhoy]) * (logT - self.temp[jT]) +
             B * (logrhoy - self.rhoy[irhoy]) + C * (logT - self.temp[jT]) + D)

        return r


class TabularRate(Rate):
    """A rate tabulated in terms of log10(ρ Y_e) and log10(T).

    Parameters
    ----------
    rfile : str, Path, io.StringIO
        the file containing the data table

    """
    def __init__(self, rfile=None):
        super().__init__()
        self.rate_eval_needs_rho = True
        self.rate_eval_needs_comp = True

        self.rfile_path = None
        self.rfile = None
        self.source = None

        if isinstance(rfile, (str, Path)):
            rfile = Path(rfile)
            self.rfile_path = _find_rate_file(rfile)
            self.source = RateSource.source(self.rfile_path.parent.name)
            self.rfile = rfile.name

        self.fname = None

        self.label = "tabular"
        self.tabular = True

        # we should initialize this somehow
        self.weak_type = ""

        if isinstance(rfile, Path):
            # read in the file, parse the different sets and store them as
            # SingleSet objects in sets[]
            f = self.rfile_path.open()
        elif isinstance(rfile, io.StringIO):
            # Set f to the io.StringIO object
            f = rfile
        else:
            f = None

        if f:
            self._read_from_file(f)
            f.close()

        self._set_rhs_properties()
        self._set_screening()
        self._set_print_representation()

        self.get_tabular_rate()

        # store the extrema of the thermodynamics
        _rhoy = self.tabular_data_table[::self.table_temp_lines, TableIndex.RHOY.value]
        _temp = self.tabular_data_table[0:self.table_temp_lines, TableIndex.T.value]

        self.table_Tmin = 10.0**(_temp.min())
        self.table_Tmax = 10.0**(_temp.max())
        self.table_rhoYmin = 10.0**(_rhoy.min())
        self.table_rhoYmax = 10.0**(_rhoy.max())

        self.interpolator = TableInterpolator(self.table_rhoy_lines, self.table_temp_lines,
                                              self.tabular_data_table)

    def __hash__(self):
        return hash(self.__repr__())

    def __eq__(self, other):
        """Determine whether two Rate objects are equal.  They are
        equal if they contain identical reactants and products.

        """

        if not isinstance(other, TabularRate):
            return False

        return self.reactants == other.reactants and self.products == other.products

    def __add__(self, other):
        raise NotImplementedError("addition not defined for tabular rates")

    def _read_from_file(self, f):
        """Given a file object, read rate data from the file.

        Parameters
        ----------
        f : io.TextIOWrapper
            The file object

        """
        lines = f.readlines()
        f.close()

        self.original_source = "".join(lines)

        # first line is the chapter
        self.chapter = lines[0].strip()
        if self.chapter != "t":
            raise RateFileError(f"Invalid chapter for TabularRate ({self.chapter})")

        # remove any blank lines
        set_lines = [line for line in lines[1:] if not line.strip() == ""]

        # e1 -> e2, Tabulated
        s1 = set_lines.pop(0)
        s2 = set_lines.pop(0)
        s3 = set_lines.pop(0)
        s4 = set_lines.pop(0)
        s5 = set_lines.pop(0)
        f = s1.split()
        try:
            self.reactants.append(Nucleus.from_cache(f[0]))
            self.products.append(Nucleus.from_cache(f[1]))
        except UnsupportedNucleus as ex:
            raise RateFileError(f'Nucleus objects could not be identified in {self.original_source}') from ex

        self.table_file = s2.strip()
        self.table_header_lines = int(s3.strip())
        self.table_rhoy_lines = int(s4.strip())
        self.table_temp_lines = int(s5.strip())
        self.table_num_vars = 6  # Hard-coded number of variables in tables for now.
        self.table_index_name = f'j_{self.reactants[0]}_{self.products[0]}'
        self.labelprops = 'tabular'

        # set weak type
        if "electroncapture" in self.table_file:
            self.weak_type = "electron_capture"

        elif "betadecay" in self.table_file:
            self.weak_type = "beta_decay"

        # since the reactants and products were only now set, we need
        # to recompute Q -- this is used for finding rate pairs
        self._set_q()

    def _set_rhs_properties(self):
        """Compute statistical prefactor and density exponent from the
        reactants."""
        self.prefactor = 1.0  # this is 1/2 for rates like a + a (double counting)
        self.inv_prefactor = 1
        if self.use_identical_particle_factor:
            for r in set(self.reactants):
                self.inv_prefactor = self.inv_prefactor * math.factorial(self.reactants.count(r))
        self.prefactor = self.prefactor/float(self.inv_prefactor)
        self.dens_exp = len(self.reactants)-1

    def _set_screening(self):
        """Tabular rates are not currently screened (they are
        e-capture or beta-decay)"""
        self.ion_screen = []
        self.symmetric_screen = []

        if not self.fname:
            # This is used to determine which rates to detect as the same reaction
            # from multiple sources in a Library file, so it should not be unique
            # to a given source, e.g. wc12, but only unique to the reaction.
            reactants_str = '_'.join([repr(nuc) for nuc in self.reactants])
            products_str = '_'.join([repr(nuc) for nuc in self.products])
            self.fname = f'{reactants_str}__{products_str}'

    def get_rate_id(self):
        """Get an identifying string for this rate.

        Returns
        -------
        str
        """

        ssrc = 'tabular'

        return f'{self.rid} <{self.label.strip()}_{ssrc}>'

    def function_string_py(self):
        """
        Construct the python function that computes the rate.

        Returns
        -------
        str

        """

        fstring = ""
        fstring += "@numba.njit()\n"
        fstring += f"def {self.fname}(rate_eval, T, rho, Y):\n"
        fstring += f"    # {self.rid}\n"
        fstring += "    rhoY = rho * ye(Y)\n"

        fstring += f"    {self.fname}_interpolator = TableInterpolator(*{self.fname}_info)\n"

        fstring += f"    r = {self.fname}_interpolator.interpolate(np.log10(rhoY), np.log10(T), TableIndex.RATE.value)\n"
        fstring += f"    rate_eval.{self.fname} = 10.0**r\n\n"

        return fstring

    def get_tabular_rate(self):
        """Read the rate data from .dat file """

        # find .dat file and read it
        self.table_path = _find_rate_file(self.table_file)
        t_data2d = []
        with self.table_path.open() as tabular_file:
            for i, line in enumerate(tabular_file):
                # skip header lines
                if i < self.table_header_lines:
                    continue
                line = line.strip()
                # skip empty lines
                if not line:
                    continue
                # split the column values on whitespace
                t_data2d.append(line.split())

        # convert the nested list of string values into a numpy float array
        self.tabular_data_table = np.array(t_data2d, dtype=np.float64)

    def eval(self, T, *, rho=None, comp=None):
        """Evaluate the reaction rate.

        Parameters
        ----------
        T : float
            the temperature to evaluate the rate at
        rho : float
            the density to evaluate the rate at (not needed for ReacLib
            rates).
        comp : float
            the composition (of type
            :py:class:`Composition <pynucastro.networks.rate_collection.Composition>`)
            to evaluate the rate with (not needed for ReacLib rates).

        Returns
        -------
        float

        """

        rhoY = rho * comp.ye
        r = self.interpolator.interpolate(np.log10(rhoY), np.log10(T),
                                          TableIndex.RATE.value)
        return 10.0**r

    def get_nu_loss(self, T, *, rho=None, comp=None):
        """Evaluate the neutrino loss for the rate.

        Parameters
        ----------
        T : float
            the temperature to evaluate the rate at
        rho : float
            the density to evaluate the rate at (not needed for ReacLib
            rates).
        comp : float
            the composition (of type
            :py:class:`Composition <pynucastro.networks.rate_collection.Composition>`)
            to evaluate the rate with (not needed for ReacLib rates).

        Returns
        -------
        float

        """

        rhoY = rho * comp.ye
        r = self.interpolator.interpolate(np.log10(rhoY), np.log10(T),
                                          TableIndex.NU.value)
        return 10**r

    def plot(self, *, Tmin=None, Tmax=None, rhoYmin=None, rhoYmax=None,
             color_field='rate', figsize=(10, 10)):
        """plot the rate's temperature sensitivity vs temperature

        :param float Tmin:    minimum temperature for plot
        :param float Tmax:    maximum temperature for plot
        :param float rhoYmin: minimum electron density to plot (e-capture rates only)
        :param float rhoYmax: maximum electron density to plot (e-capture rates only)
        :param tuple figsize: figure size specification for matplotlib

        :return: a matplotlib figure object
        :rtype: matplotlib.figure.Figure

        """

        fig, ax = plt.subplots(figsize=figsize)

        if Tmin is None:
            Tmin = self.table_Tmin
        if Tmax is None:
            Tmax = self.table_Tmax
        if rhoYmin is None:
            rhoYmin = self.table_rhoYmin
        if rhoYmax is None:
            rhoYmax = self.table_rhoYmax

        data = self.tabular_data_table

        inde1 = data[:, TableIndex.T.value] <= np.log10(Tmax)
        inde2 = data[:, TableIndex.T.value] >= np.log10(Tmin)
        inde3 = data[:, TableIndex.RHOY.value] <= np.log10(rhoYmax)
        inde4 = data[:, TableIndex.RHOY.value] >= np.log10(rhoYmin)
        data_heatmap = data[inde1 & inde2 & inde3 & inde4].copy()

        rows, row_pos = np.unique(data_heatmap[:, 0], return_inverse=True)
        cols, col_pos = np.unique(data_heatmap[:, 1], return_inverse=True)
        pivot_table = np.zeros((len(rows), len(cols)), dtype=data_heatmap.dtype)

        if color_field == 'rate':
            icol = TableIndex.RATE.value
            title = f"{self.weak_type} rate in log10(1/s)"
            cmap = 'magma'

        elif color_field == 'nu_loss':
            icol = TableIndex.NU.value
            title = "neutrino energy loss rate in log10(erg/s)"
            cmap = 'viridis'

        else:
            raise ValueError("color_field must be either 'rate' or 'nu_loss'.")

        try:
            pivot_table[row_pos, col_pos] = data_heatmap[:, icol]
        except ValueError:
            print("Divide by zero encountered in log10\nChange the scale of T or rhoY")

        im = ax.imshow(pivot_table, cmap=cmap, origin="lower",
                       extent=[np.log10(Tmin), np.log10(Tmax), np.log10(rhoYmin), np.log10(rhoYmax)])
        fig.colorbar(im, ax=ax)

        ax.set_xlabel(r"$\log(T)$ [K]")
        ax.set_ylabel(r"$\log(\rho Y_e)$ [g/cm$^3$]")
        ax.set_title(fr"{self.pretty_string}" + "\n" + title)

        return fig
