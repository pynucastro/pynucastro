"""Support modules to write a pure python reaction network ODE source."""

import io
import sys
import types
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp

from pynucastro.constants import constants
from pynucastro.networks.rate_collection import Composition, RateCollection
from pynucastro.rates import ApproximateRate, ModifiedRate
from pynucastro.screening import get_screening_func, get_screening_map


class NetworkSolution:
    """A class to hold the solution from integrating PythonNetwork.
    The member functions provide common visualization and
    analysis routines.

    Parameters
    ----------
    sol : object
        Solution object returned by :func:`scipy.integrate.solve_ivp`. The
        array `sol.y` is assumed to contain the molar abundances, `Y_i`,
        ordered consistently with `unique_nuclei`
    rhs : Callable
        Function that computes the RHS of the PythonNetwork
    jac : Callable
        Function that computes the Jacobian of the PythonNetwork
    network : PythonNetwork
        PythonNetwork used for integration
    rho : float
        density used to integrate the network
    T : float
        temperature used to integrate the network
    screen_func: Callable
        screening function used to evaluate rates when integrating
        the network

    """

    def __init__(self, sol, rhs, jac, network, rho, T, screen_func=None):

        self._sol = sol
        self._rhs = rhs
        self._jac = jac
        self.network = network
        self.rho = rho
        self.T = T
        self.screen_func = screen_func

    @property
    def t(self):
        """Return the time array for integration

        Returns
        -------
        numpy.ndarray

        """

        return self._sol.t

    @property
    def Y(self):
        """Return the 2D array of molar abundances
        for all times.

        Returns
        -------
        numpy.ndarray

        """

        return self._sol.y

    @property
    def unique_nuclei(self):
        """Return a list of nuclei explicitely carried in the network,
        ordered consistent with molar fraction solution, Y.

        Returns
        -------
        List(Nucleus)
        """

        return self.network.unique_nuclei

    def Y_at(self, t):
        """Evaluate the molar abundances for a given time.

        Parameters
        ----------
        t : float or list or numpy.ndarray
            time or time array used to evaluate the molar abundances

        Returns
        -------
        numpy.ndarray

        """

        return self._sol.sol(t)

    def ye(self, Y):
        """Evaluate the electron fraction with a given set of molar fractions

        Parameters
        ----------
        Y : numpy.ndarray
            Molar fraction array

        Returns
        -------
        float

        """

        ye = sum(nuc.Z * Y[i] for i, nuc in enumerate(self.unique_nuclei))
        return ye

    def ye_at(self, t):
        """Evaluate the electron fraction for a given time

        Parameters
        ----------
        t : float
            time used to evaluate the electron fraction

        Returns
        -------
        float

        """

        Y = self.Y_at(t)
        return self.ye(Y)

    def rhs(self, t, Y):
        """Evaluate the RHS of the network with the same thermodynamic
        condition and screening routine used to integrate the network.

        Parameters
        ----------
        t : float
            time used to evaluate the RHS
        Y : numpy.ndarray
            molar abundances of the species

        Returns
        -------
        numpy.ndarray

        """

        return self._rhs(t, Y, self.rho, self.T, screen_func=self.screen_func)

    def rhs_at(self, t):
        """Evaluate the RHS of the network for a given time.

        Parameters
        ----------
        t : float
            time used to evaluate the RHS

        Returns
        -------
        numpy.ndarray

        """

        Y = self.Y_at(t)
        return self.rhs(t, Y)

    def jac(self, t, Y):
        """Evaluate the Jacobian of the network with the same thermodynamic
        condition and screening routine used to integrate the network.

        Parameters
        ----------
        t : float
            time used to evaluate the RHS
        Y : numpy.ndarray
            molar abundances of the species

        Returns
        -------
        numpy.ndarray

        """

        return self._jac(t, Y, self.rho, self.T, screen_func=self.screen_func)

    def jac_at(self, t):
        """Evaluate the Jacobian of the network for a given time.

        Parameters
        ----------
        t : float
            time used to evaluate the Jacobian

        Returns
        -------
        numpy.ndarray

        """

        Y = self.Y_at(t)
        return self.jac(t, Y)

    def energy_release(self, dY):
        """Evaluate the energy release in erg/g (/s if dY is actually dY/dt)

        Parameters
        ----------
        dY : numpy.ndarray
            Finite change or the rate of instataneous change in molar fractions

        Returns
        -------
        float

        """

        enuc = sum(nuc.mass * dY[i] for i, nuc in enumerate(self.unique_nuclei))
        enuc *= -1*constants.N_A*constants.MeV2erg
        return enuc

    def energy_release_at(self, t):
        """Evaluate the instantaneous energy release in erg/g/s for a given time

        Parameters
        ----------
        t : float
            time used to evaluate the instantaneous energy release

        Returns
        -------
        float

        """

        dYdt = self.rhs_at(t)
        enuc = self.energy_release(dYdt)
        return enuc

    def plot_evolution(self,
                       tmin=None, tmax=None,
                       ymin=None, ymax=None,
                       size=(800, 600), dpi=100,
                       X_cutoff_value=None,
                       label_size=14, legend_size=10,
                       three_level_style=False,
                       outfile=None):
        """Plot the time evolution of nuclei mass fractions using the
        solution returned by SciPy's solve_ivp().

        Parameters
        ----------
        tmin : float
            Minimum time shown on the x-axis. If `None`, the first value of
            `self.t` is used.
        tmax : float
            Maximum time shown on the x-axis. If `None`, the last value of
            `self.t` is used.
        ymin : float
            Minimum mass fraction shown on the y-axis. If `None`,
            use the Matplotlib autoscaled value.
        ymax : float
            Maximum mass fraction shown on the y-axis. If `None`,
            use the Matplotlib autoscaled value. The autoscaled value
            is capped at 1.2
        dpi : int
            dots per inch used with size to set output image size
        size : (tuple, list)
            (width, height) of the plot in pixels
        X_cutoff_value : float
            Minimum peak mass fraction required for a nucleus to be plotted.
        label_size : int
            Font size for axis labels.
        legend_size : int
            Font size for the legend.
        three_level_style : bool
            If `True`, use three-level linestyle and linewidth based on the peak
            mass fraction to help distinguish different curves.
            If `False`, all curves use the same line style and linewidth.
        outfile : str
            output name of the plot (extension determines the type)

        Returns
        -------
        matplotlib.figure.Figure

        """

        fig, ax = plt.subplots(figsize=(size[0]/dpi, size[1]/dpi))
        for i, nuc in enumerate(self.unique_nuclei):

            X = self.Y[i, :] * nuc.A
            max_X = X.max()
            if X_cutoff_value is None and ymin is not None:
                X_cutoff_value = ymin
            if X_cutoff_value is not None and max_X <= X_cutoff_value:
                continue

            # Set linestyle and linewidth
            lw = 1.5
            ls = "-"
            if three_level_style:
                # Set 3 levels of visual levels depending on maximum mass fraction
                lw = 1
                ls = "--"
                if max_X > 0.5:
                    lw = 2.5
                    ls = "-"
                elif max_X > 0.01:
                    lw = 1.5
                    ls = "-"

            ax.loglog(self.t, X, lw=lw, ls=ls,
                      label=rf"X(${nuc.pretty}$)")

        if tmin is None:
            tmin = self.t[0]
        if tmax is None:
            tmax = self.t[-1]

        # Auto set number of legend column
        ncol = max(1, len(ax.lines) // 8 + 1)

        ax.set_xlim(tmin, tmax)
        ax.set_ylim(ymin, ymax)
        cur_ymin, cur_ymax = ax.get_ylim()
        if ymax is None and cur_ymax > 1.2:
            # Make sure the autoscaled ymax is not greater than 1.2
            cur_ymax = 1.2
        ax.set_ylim(cur_ymin, cur_ymax)

        ax.set_xlabel("time [s]", fontsize=label_size)
        ax.set_ylabel("X", fontsize=label_size)
        ax.legend(loc="best", fontsize=legend_size, ncol=ncol)
        ax.grid(ls=":")
        fig.tight_layout()

        if outfile is not None:
            fig.savefig(outfile, dpi=dpi)

        return fig


class PythonNetwork(RateCollection):
    """A pure python reaction network.  This can create a python
    module as a file that contains everything needed to evaluate the
    reaction rates and construct the righthand side and Jacobian
    functions.

    """

    def full_ydot_string(self, nucleus, indent=""):
        """Construct a string containing the python code for
        dY(nucleus)/dt by considering every reaction that involves
        nucleus, adding terms that result in its creation and
        subtracting terms representing its destruction.

        Parameters
        ----------
        nucleus : Nucleus
            The nucleus we are constructing the time derivative
            of.
        indent : str
            A string that will be prepended to each line of the output,
            typically consisting of just spaces representing the amount
            of indent desired.

        Returns
        -------
        str

        """

        ostr = ""
        if not self.nuclei_consumed[nucleus] + self.nuclei_produced[nucleus]:
            # this captures an inert nucleus
            ostr += f"{indent}dYdt[j{nucleus.raw}] = 0.0\n\n"
        else:
            ostr += f"{indent}dYdt[j{nucleus.raw}] = (\n"
            for ipair, rp in enumerate(self.nuclei_rate_pairs[nucleus]):
                # when we are working with rate pairs, one or more of the
                # rates may be missing.  We also have not clearly separated
                # them into creation / destruction, so we'll figure that out
                rlist = [r for r in [rp.forward, rp.reverse] if r is not None]
                ostr += f"{indent}      "
                if len(rlist) > 1:
                    ostr += "( "

                for rate in rlist:
                    c_reac = rate.reactant_count(nucleus)
                    c_prod = rate.product_count(nucleus)
                    c = c_prod - c_reac
                    if c == 1:
                        ostr += f"+{rate.ydot_string_py()} "
                    elif c == -1:
                        ostr += f"-{rate.ydot_string_py()} "
                    else:
                        ostr += f"+ {c}*{rate.ydot_string_py()} "

                if len(rlist) > 1:
                    ostr += ")"
                if ipair < len(self.nuclei_rate_pairs[nucleus]) - 1:
                    ostr += " +"
                ostr = ostr.rstrip() + "\n"

            ostr += f"{indent}   )\n\n"

        return ostr

    def full_jacobian_element_string(self, ydot_i_nucleus, y_j_nucleus, indent=""):
        """Construct a string containing the python code for a single
        element of the Jacobian, dYdot(ydot_i_nucleus)/dY(y_j_nucleus)

        Parameters
        ----------
        ydot_i_nucleus: Nucleus
            The nucleus representing the dY/dt term we are differentiating.
            This is the row of the Jacobian.
        ydot_j_nucleus: Nucleus
            The nucleus we are differentiating with respect to.  This
            is the column of the Jacobian.
        indent : str
            A string that will be prepended to each line of the output,
            typically consisting of just spaces representing the amount
            of indent desired.

        Returns
        -------
        str
        """

        # this is the jac(i,j) string
        idx_str = f"jac[j{ydot_i_nucleus.raw}, j{y_j_nucleus.raw}]"

        ostr = ""
        if not self.nuclei_consumed[ydot_i_nucleus] + self.nuclei_produced[ydot_i_nucleus]:
            # this covers the case where a nucleus is not created or
            # destroyed in the entire network, but is just passive
            ostr += f"{indent}{idx_str} = 0.0\n\n"
        else:
            ostr += f"{indent}{idx_str} = (\n"
            rate_terms_str = ""
            for r in self.nuclei_consumed[ydot_i_nucleus]:
                c = r.reactant_count(ydot_i_nucleus)

                jac_str = r.jacobian_string_py(y_j_nucleus)
                if jac_str == "":
                    continue

                if c == 1:
                    rate_terms_str += f"{indent}   -{jac_str}\n"
                else:
                    rate_terms_str += f"{indent}   -{c}*{jac_str}\n"
            for r in self.nuclei_produced[ydot_i_nucleus]:
                c = r.product_count(ydot_i_nucleus)

                jac_str = r.jacobian_string_py(y_j_nucleus)
                if jac_str == "":
                    continue

                if c == 1:
                    rate_terms_str += f"{indent}   +{jac_str}\n"
                else:
                    rate_terms_str += f"{indent}   +{c}*{jac_str}\n"

            if rate_terms_str == "":
                return ""
            ostr += rate_terms_str

            ostr += f"{indent}   )\n\n"

        return ostr

    def screening_string(self, indent=""):
        """Create a string containing the python code that sets up the
        screening (PlasmaState) and calls the screening function on
        every set of reactants in our network. This computes log(screening)
        terms and store them in local variables.

        Parameters
        ----------
        indent : str
            A string that will be prepended to each line of the output,
            typically consisting of just spaces representing the amount
            of indent desired.

        Returns
        -------
        str
        """

        ostr = ""

        screening_map = get_screening_map(self.get_rates())

        # Initialize log_scor to 0.0
        for scr in screening_map:
            screen_var = f"log_scor_{scr.n1}_{scr.n2}"
            ostr += f"{indent}{screen_var} = 0.0\n"

        # Check if we're doing screening, return early if not screening
        if not self.do_screening:
            return ostr

        ostr += "\n"
        ostr += f"{indent}if screen_func is not None:\n"
        indent += "    "
        ostr += f"{indent}plasma_state = PlasmaState(T, rho, Y, Z)\n\n"
        for scr in screening_map:
            screen_var = f"log_scor_{scr.n1}_{scr.n2}"
            ostr += f"{indent}scn_fac = ScreenFactors({scr.n1.Z}, {scr.n1.A}, {scr.n2.Z}, {scr.n2.A})\n"
            ostr += f"{indent}{screen_var} = screen_func(plasma_state, scn_fac)\n"

        return ostr

    def rates_string(self, indent=""):
        """Create the python code that calls the evaluation function
        for each rate.  Typically this is of the form
        ``name(rate_eval, ...)``, where ``rate_eval`` is a container
        class in the output network that stores the rate values.  This
        also calls ``screening_string()`` after the main rates are
        evaluated but before the approximate rates are constructed.

        Parameters
        ----------
        indent : str
            A string that will be prepended to each line of the output,
            typically consisting of just spaces representing the amount
            of indent desired.

        Returns
        -------
        str

        """

        def format_rate_call(r, use_tf=True):
            args = ["rate_eval"]
            if use_tf:
                args.append("tf")
            else:
                args.append("T")
            if r.rate_eval_needs_rho:
                args.append("rho=rho")
            if r.rate_eval_needs_comp:
                args.append("Y=Y")
            if r.ion_screen:
                scr_reactants = r.ion_screen.copy()
                screen_terms = []
                while len(scr_reactants) > 1:
                    a, b = scr_reactants[0], scr_reactants[1]
                    screen_terms.append(f"log_scor_{a}_{b}")

                    # merge reactants to get compound nucleus
                    scr_reactants = [a + b] + scr_reactants[2:]
                    scr_reactants.sort(key=lambda x: x.Z)

                args.append("log_scor=" + "+".join(screen_terms))
            return f"{indent}{r.fname}({', '.join(args)})\n"

        ostr = ""

        # Precompute screening terms. Note here we compute log_screening
        ostr += self.screening_string(indent=indent)
        ostr += "\n"

        ostr += f"{indent}# reaclib rates\n"
        for r in self.reaclib_rates:
            ostr += format_rate_call(r)

        if self.tabular_rates:
            ostr += f"\n{indent}# tabular rates\n"
        for r in self.tabular_rates:
            ostr += format_rate_call(r, use_tf=False)

        if self.temperature_tabular_rates:
            ostr += f"\n{indent}# temperature tabular rates\n"
        for r in self.temperature_tabular_rates:
            ostr += format_rate_call(r, use_tf=False)

        if self.custom_rates:
            ostr += f"\n{indent}# custom rates\n"
        for r in self.custom_rates:
            ostr += format_rate_call(r)

        # modified rates will have their own screening,
        # either using the original rate or any modified
        # form.  Therefore we call them before applying
        # screening factors.

        if self.modified_rates:
            ostr += f"\n{indent}# modified rates\n"
        for r in self.modified_rates:
            ostr += format_rate_call(r)

        # Derived rate should go last (before approx rates)
        # since the inverse rate should be evaluated first.
        if self.derived_rates:
            ostr += f"\n{indent}# derived rates\n"
        for r in self.derived_rates:
            ostr += format_rate_call(r)

        if self.approx_rates:
            ostr += f"\n{indent}# approximate rates\n"
        for r in self.approx_rates:
            ostr += format_rate_call(r)

        return ostr

    def _write_network(self, outfile: str | Path = None):
        """Create the entire python code representing this network.
        This includes defining the nuclei properties, evaluating
        partition functions, defining the ``RateEval`` class, defining
        the tabular rate data, writing the functions that evaluate
        each of the rates, and then calling the functions that
        construct the RHS and Jacobian.

        Parameters
        ----------
        outfile : str, Path
            The name of the output file to write to.  If this is ``None``,
            then the output is written to stdout

        """
        # pylint: disable=arguments-differ

        if outfile is None:
            of = sys.stdout
            close_file = False
        elif hasattr(outfile, "write"):
            # already a file-like object
            of = outfile
            close_file = False
        else:
            outfile = Path(outfile)
            of = outfile.open("w")
            close_file = True

        indent = 4*" "

        of.write("import numba\n")
        of.write("import numpy as np\n")
        of.write("from pynucastro.constants import constants\n")
        of.write("from numba.experimental import jitclass\n\n")

        of.write("from pynucastro.rates import (TableIndex, TableInterpolator, TabularRate,\n")
        of.write("                              TempTableInterpolator, TemperatureTabularRate,\n")
        of.write("                              Tfactors)\n")
        of.write("from pynucastro.screening import PlasmaState, ScreenFactors\n\n")

        # integer keys

        for i, n in enumerate(self.unique_nuclei):
            of.write(f"j{n.raw} = {i}\n")

        of.write(f"nnuc = {len(self.unique_nuclei)}\n\n")

        # nuclei properties

        of.write("A = np.zeros((nnuc), dtype=np.int32)\n\n")
        for n in self.unique_nuclei:
            of.write(f"A[j{n.raw}] = {n.A}\n")

        of.write("\n")

        of.write("Z = np.zeros((nnuc), dtype=np.int32)\n\n")
        for n in self.unique_nuclei:
            of.write(f"Z[j{n.raw}] = {n.Z}\n")

        # we'll compute the masses here in erg

        of.write("\n")

        of.write("# masses in ergs\n")
        of.write("mass = np.zeros((nnuc), dtype=np.float64)\n\n")
        for n in self.unique_nuclei:
            mass = n.A_nuc * constants.m_u_MeV_C18 * constants.MeV2erg
            of.write(f"mass[j{n.raw}] = {mass}\n")

        of.write("\n")

        of.write("names = []\n")
        for n in self.unique_nuclei:
            name = n.short_spec_name
            if name != "n":
                name = name.capitalize()
            of.write(f"names.append(\"{name}\")\n")

        of.write("\n")

        of.write("def to_composition(Y):\n")
        of.write(f'{indent}''"""Convert an array of molar fractions to a Composition object."""\n')
        of.write(f'{indent}'"from pynucastro import Composition, Nucleus\n")
        of.write(f'{indent}'"nuclei = [Nucleus.from_cache(name) for name in names]\n")
        of.write(f'{indent}'"comp = Composition(nuclei)\n")
        of.write(f'{indent}'"for i, nuc in enumerate(nuclei):\n")
        of.write(f'{indent*2}'"comp.X[nuc] = Y[i] * A[i]\n")
        of.write(f'{indent}'"return comp\n\n")

        of.write("\n")

        of.write("def energy_release(dY):\n")
        of.write(f'{indent}''"""return the energy release in erg/g (/s if dY is actually dY/dt)"""\n')
        of.write(f'{indent}'"enuc = 0.0\n")
        of.write(f'{indent}'"for i, y in enumerate(dY):\n")
        of.write(f'{indent*2}'"enuc += y * mass[i]\n")
        of.write(f'{indent}'"enuc *= -1*constants.N_A\n")
        of.write(f'{indent}'"return enuc\n\n")

        # partition function data (if needed)

        nuclei_pfs = self.get_nuclei_needing_partition_functions()

        for n in nuclei_pfs:
            of.write(f"{n}_temp_array = np.array({list(n.partition_function.T9_points)})\n")
            of.write(f"{n}_log_pf_array = np.array({list(n.partition_function.log_pf_data)})\n")
            of.write("\n")

        # rate_eval class

        of.write("@jitclass([\n")
        for r in self.all_rates:
            of.write(f'{indent}("{r.fname}", numba.float64),\n')
        of.write("])\n")
        of.write("class RateEval:\n")
        of.write(f"{indent}def __init__(self):\n")
        for r in self.all_rates:
            of.write(f"{indent*2}self.{r.fname} = np.nan\n")

        of.write("\n")

        # tabular rate data
        if self.tabular_rates:
            of.write("# note: we cannot make the TableInterpolator global, since numba doesn't like global jitclass\n")

        for r in self.tabular_rates:
            of.write(f"# load data for {r.rid}\n")
            of.write(f"{r.fname}_info = (\n")
            of.write(f"    {r.table_rhoy_lines},    # table_rhoy_lines\n")
            of.write(f"    {r.table_temp_lines},    # table_temp_lines\n")
            of.write("    # tabular_data_table\n")
            of.write(f"    np.array({r.tabular_data_table.tolist()})\n")
            of.write(")\n\n")

        # temperature tabular rate data
        if self.temperature_tabular_rates:
            of.write("# note: we cannot make the TempTableInterpolator global, since numba doesn't like global jitclass\n")

        for r in self.temperature_tabular_rates:

            of.write(f"# temperature / rate tabulation for {r.rid}\n")

            log_temp_str = np.array2string(r.log_t9_data,
                                           max_line_width=70, precision=17, separator=", ")
            of.write(f"{r.fname}_log_t9_data = np.array(\n")
            for line in log_temp_str.split("\n"):
                of.write(f"     {line}\n")
            of.write("   )\n")

            log_rate_str = np.array2string(r.log_rate_data,
                                           max_line_width=70, precision=17, separator=", ")
            of.write(f"{r.fname}_log_rate_data = np.array(\n")
            for line in log_rate_str.split("\n"):
                of.write(f"     {line}\n")
            of.write("   )\n")

            of.write(f"{r.fname}_info = ({r.fname}_log_t9_data, {r.fname}_log_rate_data)\n\n")

        # Ye helper function
        of.write("@numba.njit()\n")
        of.write("def ye(Y):\n")
        of.write(f"{indent}return np.sum(Z * Y)/np.sum(A * Y)\n\n")

        # the functions to evaluate the temperature dependence of the rates

        _rate_func_written = []
        for r in self.rates:
            if isinstance(r, ApproximateRate):
                # write out the function string for all of the rates we depend on
                for cr in r.get_child_rates():
                    if cr in _rate_func_written:
                        continue
                    of.write(cr.function_string_py())
                    _rate_func_written.append(cr)

                # now write out the function that computes the
                # approximate rate
                of.write(r.function_string_py())
            elif isinstance(r, ModifiedRate):
                orig_rate = r.original_rate
                if r in _rate_func_written:
                    continue
                of.write(orig_rate.function_string_py())
                _rate_func_written.append(orig_rate)

                # now write out the function that computes the
                # modified rate
                of.write(r.function_string_py())
                _rate_func_written.append(r)
            else:
                if r in _rate_func_written:
                    continue
                of.write(r.function_string_py())
                _rate_func_written.append(r)

        # the rhs() function

        of.write("def rhs(t, Y, rho, T, screen_func=None):\n")
        of.write(f"{indent}return rhs_eq(t, Y, rho, T, screen_func)\n\n")

        of.write("@numba.njit()\n")
        of.write("def rhs_eq(t, Y, rho, T, screen_func):\n\n")

        # get the rates
        of.write(f"{indent}tf = Tfactors(T)\n")
        of.write(f"{indent}rate_eval = RateEval()\n\n")

        of.write(self.rates_string(indent=indent))

        of.write("\n")

        of.write(f"{indent}dYdt = np.zeros((nnuc), dtype=np.float64)\n\n")

        # now make the RHSs
        for n in self.unique_nuclei:
            of.write(self.full_ydot_string(n, indent=indent))

        of.write(f"{indent}return dYdt\n\n")

        # the jacobian() function

        of.write("def jacobian(t, Y, rho, T, screen_func=None):\n")
        of.write(f"{indent}return jacobian_eq(t, Y, rho, T, screen_func)\n\n")

        of.write("@numba.njit()\n")
        of.write("def jacobian_eq(t, Y, rho, T, screen_func):\n\n")

        # get the rates
        of.write(f"{indent}tf = Tfactors(T)\n")
        of.write(f"{indent}rate_eval = RateEval()\n\n")

        of.write(self.rates_string(indent=indent))

        of.write("\n")

        of.write(f"{indent}jac = np.zeros((nnuc, nnuc), dtype=np.float64)\n\n")

        # now fill each Jacobian element
        for n_i in self.unique_nuclei:
            for n_j in self.unique_nuclei:
                of.write(self.full_jacobian_element_string(n_i, n_j, indent=indent))

        of.write(f"{indent}return jac\n")

        if close_file:
            of.close()

    def integrate_network(self, tmax, rho, T, Y0=None,
                          screen_method=None,
                          initial_comp="uniform",
                          rtol=1e-8, atol=1e-8):
        """Integrate the network to tmax given (rho, T, Y0) using
        SciPy's solve_ivp() with BDF method.

        Parameters
        ----------
        tmax: float
            final integration time.
        rho : float
            density used to integrate the network
        T : float
            temperature used to integrate the network
        Y0 : numpy.ndarray
            initial molar abundance of the nuclei. If not provided,
            the initial composition is initialized according to `initial_comp`
        screen_method : str
            name of the screening function used to evaluate rates when integrating
            the network. Valid choices are: `screen5`, `chugunov_2007`, `chugunov_2009`,
            `potekhin_1998`, and `debye_huckel`. If `None`, no screening is applied.
        initial_comp : str
            different modes to use to set up the initial composition if Y0 is None.
            Valid choices are: `uniform`, `random`, and `solar`.
        rtol : float
            relative tolerance for SciPy's solve_ivp()
        atol : float
            absolute tolerance for SciPy's solve_ivp()

        Returns
        -------
        sol : object
            Solution returned by :func:`scipy.integrate.solve_ivp`.

        """

        # Write the network module as a string
        f = io.StringIO()
        self.write_network(outfile=f)
        network_code = f.getvalue()

        # Create a new in-memory module called `network`
        network = types.ModuleType("network")

        # Execute the code inside the module namespace
        exec(network_code, network.__dict__)  # pylint: disable=exec-used

        # Get RHS and Jacobian. Use getattr to avoid pylint warning.
        rhs = getattr(network, "rhs")
        jacobian = getattr(network, "jacobian")

        # Get the appropriate screening function
        screen_func = get_screening_func(screen_method)

        # Setup the initial molar abundance if Y0 is None
        if Y0 is None:
            if initial_comp is None:
                raise ValueError("Valid initial compositions are ['uniform', 'random', 'solar']")
            comp = Composition(self.unique_nuclei, init=initial_comp)
            ys = comp.get_molar()
            Y0 = np.array([ys[nuc] for nuc in self.unique_nuclei])

        # Integrate using SciPy's solve_ivp() using BDF method -- good for stiff system.
        sol = solve_ivp(rhs, [0, tmax], Y0, method="BDF",
                        dense_output=True, args=(rho, T, screen_func),
                        rtol=rtol, atol=atol, jac=jacobian)

        # Create NetworkSolution
        network_sol = NetworkSolution(sol, rhs, jacobian, self,
                                      rho, T, screen_func=screen_func)

        return network_sol
