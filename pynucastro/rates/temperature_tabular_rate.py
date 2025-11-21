"""Classes and methods for describing a reaction rate that is
tabulated in terms of electron density and temperature.

"""

import matplotlib.pyplot as plt
import numpy as np

import pynucastro.numba_util as numba
from pynucastro.numba_util import jitclass
from pynucastro.rates import Rate


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
        return max(0, min(max_idx, np.searchsorted(self.temp_points, T9_0)) - 1)

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

        r = self.interpolator.interpolate(T)

        scor = 1.0
        if screen_func is not None:
            if rho is None or comp is None:
                raise ValueError("rho (density) and comp (Composition) needs to be defined when applying electron screening.")
            scor = self.evaluate_screening(rho, T, comp, screen_func)

        r *= scor

        return r

    def plot(self, *, Tmin=None, Tmax=None, figsize=(6, 6),
             rho=None, comp=None, screen_func=None):
        """Plot the rate as a function of temperature.

        Parameters
        ----------
        Tmin : float
            minimum temperature for the plot
        Tmax : float
            maximum temperature for the plot
        figsize : tuple
            the horizontal, vertical size (in inches) for the plot
        rho : float
            the density to evaluate the screening effect.
        comp : float
            the composition (of type
            :py:class:`Composition <pynucastro.networks.rate_collection.Composition>`)
            to evaluate the screening effect.
        screen_func : Callable
            one of the screening functions from :py:mod:`pynucastro.screening`
            -- if provided, then the rate will include the screening correction.

        Returns
        -------
        matplotlib.figure.Figure

        """

        fig, ax = plt.subplots(figsize=figsize)

        if Tmin is None:
            Tmin = self.table_Tmin
        if Tmax is None:
            Tmax = self.table_Tmax

        # pylint: disable=duplicate-code
        temps = np.logspace(np.log10(Tmin), np.log10(Tmax), 100)
        r = np.zeros_like(temps)

        for n, T in enumerate(temps):
            r[n] = self.eval(T, rho=rho, comp=comp, screen_func=screen_func)

        ax.loglog(temps, r)
        ax.set_xlabel(r"$T$")

        if self.dens_exp == 0:
            ax.set_ylabel(r"$\tau$")
        elif self.dens_exp == 1:
            ax.set_ylabel(r"$N_A <\sigma v>$")
        elif self.dens_exp == 2:
            ax.set_ylabel(r"$N_A^2 <n_a n_b n_c v>$")

        ax.set_title(fr"{self.pretty_string}")
        #pylint: enable=duplicate-code
        ax.grid(ls=":")

        return fig
