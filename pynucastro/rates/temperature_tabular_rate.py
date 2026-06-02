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
    methods that allow us to interpolate the rate.
    This uses the monotone cubic hermite interpolation for
    log(rate) and log(T9).

    Parameters
    ----------
    log_temp_points : numpy.ndarray
        an array giving the temperature at the points where we
        tabulate the rate --- this is log(T9)
    log_rate_data : numpy.ndarray
        an array giving the tabulated log(rate) data

    """

    def __init__(self, log_temp_points, log_rate_data):

        self.log_temp_points = log_temp_points
        self.log_rate_data = log_rate_data

    def _get_logT9_idx(self, log_T9_0):
        """Find the index into the temperatures such that T[i] < T0
        <= T[i+1].  We return i here, corresponding to the lower
        value. We also make sure that i and i+1 are in bounds.

        Parameters
        ----------
        log_T9_0 : float
            temperature (log(T/1.e9 K)) to interpolate at

        Returns
        -------
        int

        """

        max_idx = len(self.log_temp_points) - 2
        return max(0, min(max_idx, np.searchsorted(self.log_temp_points, log_T9_0) - 1))

    def _limit_slope(self, hm, hp, dm, dp):
        """Given the slope from the left, dm, and the slope
        from the right, dp. Limit the slope to preserve monotonicity.
        Uses harmonic mean of the two slopes weighted by the
        lengths of the two intervals.

        See C. Moler, Numerical Computing with Matlab, 2004.
        :doi:`10.1137/1.9780898717952`, Chapter 3.4

        Parameters
        ----------
        hm : float
            length interval between i and i-1
        hp : float
            length interval between i and i+1
        dm : float
            slope using i and i-1 data
        dp : float
            slope using i and i+1 data

        Returns
        -------
        float

        """

        if dm * dp <= 0.0:
            return 0.0

        w1 = 2.0 * hp + hm
        w2 = hp + 2.0 * hm
        return (w1 + w2) / (w1 / dm + w2 / dp)

    def _linear_extrap(self, T, idx0, idx1):
        """Do linear extrapolation given the log(T9).
        This returns log(rate).

        Parameters
        ----------
        T : float
            temperature (log(T/1.e9 K)) to interpolate at
        idx0 : int
            reference point index
        idx1 : int
            neighboring index for computing derivative

        Returns
        -------
        float

        """

        x0 = self.log_temp_points[idx0]
        x1 = self.log_temp_points[idx1]

        f0 = self.log_rate_data[idx0]
        f1 = self.log_rate_data[idx1]

        dfdx = (f1 - f0) / (x1 - x0)
        return f0 + dfdx * (T - x0)

    def _monotone_cubic_interp(self, T, idx):
        """Given log(T9) and the corresponding index where
        T[i] < T <= T[i+1], do Monotone Hermite cubic interpolation.
        Note this method only guarantees monotonicity within
        the table interval. It returns log(rate).

        See C. Moler, Numerical Computing with Matlab, 2004.
        :doi:`10.1137/1.9780898717952`, Chapter 3.6

        Parameters
        ----------
        T : float
            temperature (log(T/1.e9 K)) to interpolate at
        idx : int
            temperature index where T[i] < T <= T[i+1]

        Returns
        -------
        float

        """

        max_idx = len(self.log_temp_points) - 1

        # Get surround 4 data points
        # Note idx should already be clamped within [0, max_idx-1]
        assert idx >= 0 and idx < max_idx

        im1 = max(idx - 1, 0)
        ip1 = idx + 1
        ip2 = min(idx + 2, max_idx)

        x_im1 = self.log_temp_points[im1]
        x_idx = self.log_temp_points[idx]
        x_ip1 = self.log_temp_points[ip1]
        x_ip2 = self.log_temp_points[ip2]

        f_im1 = self.log_rate_data[im1]
        f_idx = self.log_rate_data[idx]
        f_ip1 = self.log_rate_data[ip1]
        f_ip2 = self.log_rate_data[ip2]

        # Compute the slope at idx and idx + 1
        m0 = m1 = (f_ip1 - f_idx) / (x_ip1 - x_idx)

        # If NOT on the edge, compute the slope from the other side
        # and limit the slope to preserve monotonicity
        if idx != 0:
            d_im1 = (f_idx - f_im1) / (x_idx - x_im1)
            m0 = self._limit_slope(x_idx - x_im1, x_ip1 - x_idx,
                                   d_im1, m0)
        if idx + 1 != max_idx:
            d_ip1 = (f_ip2 - f_ip1) / (x_ip2 - x_ip1)
            m1 = self._limit_slope(x_ip1 - x_idx, x_ip2 - x_ip1,
                                   m1, d_ip1)

        # Compute cubic Hermite polynomial basis
        h = x_ip1 - x_idx
        t = (T - x_idx) / h

        H0 = 2.0*t**3 - 3.0*t**2 + 1.0
        H1 = -2.0*t**3 + 3.0*t**2
        Hhat0 = h * (t**3 - 2.0*t**2 + t)
        Hhat1 = h * (t**3 - t**2)

        # Construct the interpolant
        # H_3(x) = f0 H_{1,0}(x) + f1 H_{1,1}(x) + f'0 Ĥ_{1,0}(x) + f'1 Ĥ_{1,0}(x)

        return f_idx * H0 + f_ip1 * H1 + m0 * Hhat0 + m1 * Hhat1

    def interpolate(self, T0):
        """Given T0, the temperature where we want the rate,
        do monotone Hermite cubic interpolation to find
        the value of the rate in the table. Use linear extrapolation
        when T0 is out of bound of the table.
        Note this returns log(rate).

        Parameters
        ----------
        T0 : float
            temperature to interpolate at

        Returns
        -------
        float

        """

        log_T9_0 = np.log(T0 * 1.e-9)
        max_idx = len(self.log_temp_points) - 1

        # Get temperature index such that T[i] < T0 <= T[i+1]
        # We ensure i is within [0, max_idx - 1]
        idx = self._get_logT9_idx(log_T9_0)

        # For temperature out of the bound, use linear extrapolation
        if log_T9_0 < self.log_temp_points[0]:
            return self._linear_extrap(log_T9_0, 0, 1)

        if log_T9_0 > self.log_temp_points[-1]:
            return self._linear_extrap(log_T9_0, max_idx, max_idx-1)

        # Use monotone Hermite cubic interpolation
        return self._monotone_cubic_interp(log_T9_0, idx)


class TemperatureTabularRate(Rate):
    """A rate whose temperature dependence is tabulated.

    Parameters
    ----------
    log_t9_data : numpy.ndarray
        The temperature (in log(T / 1.e9 K)) where we tabulate the rate
    log_rate_data : numpy.ndarray
        The tabulated log(rate) data, N_A <σv>

    """

    def __init__(self, log_t9_data, log_rate_data, rate_source=None,
                 label="temptab", **kwargs):
        super().__init__(label=label, rate_source=rate_source, **kwargs)

        self.tabular = True

        self.log_t9_data = log_t9_data
        self.log_rate_data = log_rate_data

        # T9 should be monotonically increasing
        assert np.all(log_t9_data[1:] > log_t9_data[:-1])

        # same number of data points in T and rate
        assert len(self.log_t9_data) == len(self.log_rate_data)

        # store the extrema of the thermodynamics
        self.table_Tmin = 1.e9 * np.exp(self.log_t9_data.min())
        self.table_Tmax = 1.e9 * np.exp(self.log_t9_data.max())

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

    def function_string_py(self):
        """Construct the python function that computes the rate.

        Returns
        -------
        str

        """

        fstring = ""
        fstring += "@numba.njit()\n"
        fstring += f"def {self.fname}(rate_eval, T, log_scor=0.0):\n"
        fstring += f"    # {self.rid}\n"
        fstring += f"    {self.fname}_interpolator = TempTableInterpolator(*{self.fname}_info)\n"
        fstring += f"    log_r = {self.fname}_interpolator.interpolate(T)\n"
        fstring += f"    rate_eval.{self.fname} = np.exp(log_r + log_scor)\n\n"

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

        # pylint: disable=duplicate-code
        if extra_args is None:
            extra_args = ()

        args = ["const tf_t& tfactors",
                f"const {dtype} log_scor", f"const {dtype} dlog_scor_dT",
                f"{dtype}& rate", f"{dtype}& drate_dT", *extra_args]
        fstring = ""
        fstring += "template <int do_T_derivatives>\n"
        fstring += f"{specifiers}\n"
        fstring += f"void rate_{self.fname}({', '.join(args)}) {{\n\n"
        fstring += f"    // {self.rid}\n\n"
        # pylint: enable=duplicate-code

        fstring += "    auto [_log_rate, _dlog_rate_dlogT9] = interp_net::monotone_1d_interp<do_T_derivatives>(\n"
        fstring += "                                               tfactors.lnT9,\n"
        fstring += f"                                               {self.fname}_data::log_t9,\n"
        fstring += f"                                               {self.fname}_data::log_rate);\n"
        fstring += "    rate = std::exp(_log_rate + log_scor);\n"
        fstring += "    // we found dlog(rate)/dlog(T9)\n"
        fstring += "    if constexpr (do_T_derivatives) {\n"
        fstring += f"        {dtype} dlog_rate_dT = tfactors.T9i * _dlog_rate_dlogT9 * 1.0e-9_rt + dlog_scor_dT\n;"
        fstring += "        drate_dT = rate * dlog_rate_dT;\n"
        fstring += "    }\n"

        if not leave_open:
            fstring += "}\n\n"

        return fstring

    def log_eval(self, T, *, rho=None, comp=None,
                 screen_func=None):
        """Evaluate the natural log of reaction rate for temperature T.

        Parameters
        ----------
        T : float
            the temperature to evaluate the rate at
        rho : float
            the density to evaluate the rate at.
        comp : float
            the composition (of type
            :py:class:`Composition <pynucastro.nucdata.composition.Composition>`)
            to evaluate the rate with.
        screen_func : Callable
            one of the screening functions from :py:mod:`pynucastro.screening`
            -- if provided, then the rate will include screening correction.

        Returns
        -------
        float

        """

        log_rate = self.interpolator.interpolate(T)

        log_scor = 0.0
        if screen_func is not None:
            if rho is None or comp is None:
                raise ValueError("rho (density) and comp (Composition) needs to be defined when applying electron screening.")
            log_scor = self.evaluate_screening(rho, T, comp, screen_func)

        log_rate += log_scor

        return log_rate

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
            :py:class:`Composition <pynucastro.nucdata.composition.Composition>`)
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

        temps = np.logspace(np.log10(Tmin), np.log10(Tmax), 100)
        r = np.zeros_like(temps)

        for n, T in enumerate(temps):
            r[n] = self.eval(T, rho=rho, comp=comp, screen_func=screen_func)

        ax.loglog(temps, r)

        # overplot tabular data within plotting range
        T_data = np.exp(self.log_t9_data) * 1e9
        rate_data = np.exp(self.log_rate_data)
        mask = (T_data >= Tmin) & (T_data <= Tmax)

        ax.plot(T_data[mask], rate_data[mask],
                marker='^', markersize=6, color='k',
                linestyle="none", label="Tabular Data")

        ax.set_xlabel(r"$T [K]$")
        if self.dens_exp == 0:
            ax.set_ylabel(r"$\tau$")
        elif self.dens_exp == 1:
            ax.set_ylabel(r"$N_A \langle \sigma v\rangle$")
        elif self.dens_exp == 2:
            ax.set_ylabel(r"$N_A^2 \langle n_a n_b n_c v\rangle$")

        ax.set_title(fr"{self.pretty_string}")
        ax.grid(ls=":")

        return fig
