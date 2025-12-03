"""Classes and methods for describing a reaction rate that is
tabulated in terms of electron density and temperature.

"""

import matplotlib.pyplot as plt
import numpy as np

import pynucastro.numba_util as numba
from pynucastro.numba_util import jitclass
from pynucastro.rates.rate import Rate


@jitclass([
    ('log_temp_points', numba.float64[:]),
    ('log_rate_data', numba.float64[:])
])
class TempTableInterpolator:
    """A class that holds a pointer to the rate data and
    methods that allow us to interpolate the rate

    Parameters
    ----------
    log_temp_points : numpy.ndarray
        an array giving the temperature at the points where we
        tabulate the rate --- this is log10(T9)
    log_rate_data : numpy.ndarray
        an array giving the tabulated log10(rate) data

    """

    def __init__(self, log_temp_points, log_rate_data):

        self.log_temp_points = log_temp_points
        self.log_rate_data = log_rate_data

    def _get_logT9_idx(self, log_T9_0):
        """Find the index into the temperatures such that T[i-1] < T0
        <= T[i].  We return i-1 here, corresponding to the lower
        value.  We also make sure that i-2 and i+1 are in bounds.

        Parameters
        ----------
        log_T9_0 : float
            temperature (log10(T/1.e9 K)) to interpolate at

        Returns
        -------
        int

        """

        max_idx = len(self.log_temp_points) - 2
        return max(1, min(max_idx, np.searchsorted(self.log_temp_points, log_T9_0)) - 1)

    def interpolate(self, T0):
        """Given T0, the temperature where we want the rate, do
        cubic interpolation to find the value of the rate in the
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
        log_T9_0 = np.log10(T9_0)

        # we'll give a little epsilon buffer here to allow for roundoff
        eps = 0.005
        if log_T9_0 < self.log_temp_points.min() - eps or log_T9_0 > self.log_temp_points.max() + eps:
            raise ValueError("temperature out of table bounds")

        idx_t = self._get_logT9_idx(log_T9_0)

        # get the 4 points surrounding T0

        xp = np.array([self.log_temp_points[idx_t-1],
                       self.log_temp_points[idx_t],
                       self.log_temp_points[idx_t+1],
                       self.log_temp_points[idx_t+2]])

        fp = np.array([self.log_rate_data[idx_t-1],
                       self.log_rate_data[idx_t],
                       self.log_rate_data[idx_t+1],
                       self.log_rate_data[idx_t+2]])

        r = 0
        for m in range(4):

            # create the Lagrange basis function for point m
            l = 1
            for n in range(4):
                if n == m:
                    continue

                l *= (log_T9_0 - xp[n]) / (xp[m] - xp[n])

            r += fp[m] * l

        return 10.0**r


class TemperatureTabularRate(Rate):
    """A rate whose temperature dependence is tabulated.

    Note: presently this only supports strong-mediated rates.

    Parameters
    ----------
    log_t9_data : numpy.ndarray
        The temperature (in log10(T / 1.e9 K)) where we tabulate the rate
    log_rate_data : numpy.ndarray
        The tabulated log10(rate) data, N_A <σv>

    """

    def __init__(self, log_t9_data, log_rate_data, **kwargs):
        super().__init__(**kwargs)

        # make sure there are no weak interactions -- we don't
        # support those yet
        assert sum(n.Z for n in self.reactants) == sum(n.Z for n in self.products)
        assert sum(n.A for n in self.reactants) == sum(n.A for n in self.products)

        self.log_t9_data = log_t9_data
        self.log_rate_data = log_rate_data

        # T9 should be monotonically increasing
        assert np.all(log_t9_data[1:] > log_t9_data[:-1])

        # same number of data points in T and rate
        assert len(self.log_t9_data) == len(self.log_rate_data)

        reactants_str = '_'.join([repr(nuc) for nuc in self.reactants])
        products_str = '_'.join([repr(nuc) for nuc in self.products])
        self.fname = f'{reactants_str}__{products_str}__temptab'

        self.label = "temptab"

        # tabular here really means weak rate tabular
        self.tabular = False

        # we should initialize this somehow
        self.weak_type = ""

        # this should not really be needed
        self.chapter = "temp_tabular"

        self._set_rhs_properties()
        self._set_screening()
        self._set_print_representation()

        # store the extrema of the thermodynamics
        self.table_Tmin = 1.e9 * 10.0**self.log_t9_data.min()
        self.table_Tmax = 1.e9 * 10.0**self.log_t9_data.max()

        self.interpolator = TempTableInterpolator(self.log_t9_data, self.log_rate_data)

    def __eq__(self, other):
        """Determine whether two Rate objects are equal.  They are
        equal if they contain identical reactants and products.

        """

        if not isinstance(other, TemperatureTabularRate):
            return False

        return self.reactants == other.reactants and self.products == other.products

    def __hash__(self):
        return hash(self.__repr__())

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

        fstring = ""
        fstring += "@numba.njit()\n"
        fstring += f"def {self.fname}(rate_eval, T):\n"
        fstring += f"    # {self.rid}\n"
        fstring += f"    {self.fname}_interpolator = TempTableInterpolator(*{self.fname}_info)\n"

        fstring += f"    r = {self.fname}_interpolator.interpolate(T)\n"
        fstring += f"    rate_eval.{self.fname} = r\n\n"

        return fstring

    def function_string_cxx(self, dtype="double", specifiers="inline",
                            leave_open=False, extra_args=None):
        """Return a string containing the C++ function that computes
        the rate

        Parameters
        ----------
        dtype : str
            The C++ datatype to use for all declarations
        specifiers : str
            C++ specifiers to add before each function declaration
            (i.e. "inline")
        leave_open : bool
            If ``true``, then we leave the function unclosed (no "}"
            at the end).  This can allow additional functions to add
            to this output.
        extra_args : list, tuple
            A list of strings representing additional arguments that
            should be appended to the argument list when defining the
            function interface.

        Returns
        -------
        str

        """

        if extra_args is None:
            extra_args = ()

        args = ["const tf_t& tfactors", f"{dtype}& rate", f"{dtype}& drate_dT", *extra_args]
        fstring = ""
        fstring += "template <int do_T_derivatives>\n"
        fstring += f"{specifiers}\n"
        fstring += f"void rate_{self.cname()}({', '.join(args)}) {{\n\n"
        fstring += f"    // {self.rid}\n\n"
        fstring += "    amrex::Real log_t9 = tfactors.lnT9 * ln10_inv;\n"
        fstring += "    auto [_rate, _drate_dT] = interp_net::cubic_interp_uneven<do_T_derivatives>(\n"
        fstring += "                                               log_t9,\n"
        fstring += f"                                               {self.cname()}_data::log_t9,\n"
        fstring += f"                                               {self.cname()}_data::log_rate);\n"
        fstring += "    rate = std::pow(10.0_rt, _rate);\n"
        fstring += "    // we found dlog10(rate)/dlog10(T9)\n"
        fstring += "    if constexpr (do_T_derivatives) {\n"
        fstring += "        drate_dT = (rate / tfactors.T9) * _drate_dT * 1.e-9;\n"
        fstring += "    }\n"

        if not leave_open:
            fstring += "}\n\n"

        return fstring

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
