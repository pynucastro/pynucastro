"""Classes and methods for describing a reaction rate that is
tabulated in terms of electron density and temperature.

"""

import math
import re
from enum import Enum
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

import pynucastro.numba_util as numba
from pynucastro.nucdata import Nucleus, UnsupportedNucleus
from pynucastro.numba_util import jitclass
from pynucastro.rates.files import RateFileError, _find_rate_file
from pynucastro.rates.rate import Rate, RateSource



@jitclass([
    ('temp_points', numba.float64[:]),
    ('rate_data', numba.float64[:])
])
class TableInterpolator:
    """A class that holds a pointer to the rate data and
    methods that allow us to interpolate the rate

    Parameters
    ----------
    temp_points : numpy.ndarray
        an array giving the temperature at the points where we
        tabulate the rate --- this is assumed to be T9
    rate_data : numpy.ndarray
        an array giving the tabulated rate data

    """

    def __init__(self, temp_points, rate_data):

        self.temp_points = temp_points
        self.rate_data = rate_data

    def _get_T_idx(self, T0):
        """Find the index into the temperatures such that T[i-1] < T0 <= T[i].
        We return i-1 here, corresponding to the lower value.

        Parameters
        ----------
        T0 : float
            temperature to interpolate at

        Returns
        -------
        int

        """

        T9_0 = T0 / 1.e9

        max_idx = len(self.temp_points) - 1
        return max(0, min(max_idx, np.searchsorted(self.temp, T9_0)) - 1)

    def interpolate(self, T0):
        """Given T0, the temperature where we want the rate, do
        linear interpolation to find the value of the rate in the
        table

        Parameters
        ----------
        T0 : float
            temperature to interpolate at

        Returns
        -------
        float

        """

        T9_0 = T0 / 1.e9

        if T9_0 < self.temp_points.min() or T9_0 > self.temp_points.max():
            raise ValueError("temperature out of table bounds")


        idx_t = self._get_T_idx(T0)

        # note: T is stored as T9

        dT9 = self.temp_points[idx_t+1] - self.temp_points[idx_t]

        rate_i = self.rate_data[idx_t]
        rate_ip1 = self.rate_data[idx_t+1]

        r = rate_i + (rate_ip1 - rate_i) / dT9 * (T9_0 - self.temp_points[idx_t])
        return r


class TemperatureTabularRate(Rate):
    """A rate whose temperature dependence is tabulated.

    Note: presently this only supports strong-mediated rates.

    Parameters
    ----------
    t9_data : numpy.ndarray
        The temperature (in 1.e9 K) where we tabulate the rate
    rate_data : numpy.ndarray
        The tabulated rate data, N_A <σv>

    """

    def __init__(self, t9_data, rate_data):
        super().__init__()

        # make sure there are no weak interactions -- we don't
        # support those yet
        assert sum(n.Z for n in self.reactants) == sum(n.Z for n in self.products)
        assert sum(n.A for n in self.reactants) == sum(n.A for n in self.products)

        self.t9_data = t9_data
        self.rate_data = rate_data

        self.fname = None

        self.label = "temp_tabular"

        # tabular here really means weak rate tabular
        self.tabular = False

        # we should initialize this somehow
        self.weak_type = ""


        self._set_rhs_properties()
        self._set_screening()
        self._set_print_representation()

        # store the extrema of the thermodynamics
        self.table_Tmin = 10.0**(self.t9_data.min())
        self.table_Tmax = 10.0**(self.t9_data.max())

        self.interpolator = TableInterpolator(self.t9_data, self.rate_data)

    def __eq__(self, other):
        """Determine whether two Rate objects are equal.  They are
        equal if they contain identical reactants and products.

        """

        if not isinstance(other, TemperatureTabularRate):
            return False

        return self.reactants == other.reactants and self.products == other.products

    def __add__(self, other):
        raise NotImplementedError("addition not defined for tabular rates")

    def get_rate_id(self):
        """Get an identifying string for this rate.

        Returns
        -------
        str

        """

        ssrc = 'temp_tabular'

        return f'{self.rid} <{self.label.strip()}_{ssrc}>'

    def function_string_py(self):
        """Construct the python function that computes the rate.

        Returns
        -------
        str

        """

        raise NotImplementedError()

    def eval(self, T, *, rho=None, comp=None,
             screen_func=None):
        """Evaluate the reaction rate.

        Parameters
        ----------
        T : float
            the temperature to evaluate the rate at
        rho : float
            the density to evaluate the rate at.
        comp : float
            the composition (of type
            :py:class:`Composition <pynucastro.networks.rate_collection.Composition>`)
            to evaluate the rate with.
        screen_func : Callable
            one of the screening functions from :py:mod:`pynucastro.screening`
            -- if provided, then the rate will include screening correction.
            Unused for Tabular weak rates since screening does not affect weak reactions.

        Returns
        -------
        float

        """

        r = self.interpolator.interpolate(np.log10(rhoY), np.log10(T),
                                          TableIndex.RATE.value)
        return 10.0**r

    def plot(self, *, Tmin=None, Tmax=None, rhoYmin=None, rhoYmax=None,
             color_field='rate', figsize=(10, 10)):
        """Plot the rate or neutrino loss in the log10(ρ Y_e) and
        log10(T) plane.

        Parameters
        ----------
        Tmin : float
            minimum temperature for the plot
        Tmax : float
            maximum temperature for the plot
        rhoYmin : float
            minimum (ρ Y_e) for the plto
        rhoYmax : float
            maximum (ρ Y_e) for the plto
        color_field : str
            the field to plot.  Possible values are "rate" or "nu_loss"
        figsize : tuple
            the horizontal, vertical size (in inches) for the plot

        Returns
        -------
        matplotlib.figure.Figure

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
