"""Support modules to write a pure python reaction network ODE
source"""

import shutil
import sys
import warnings
from pathlib import Path

from pynucastro.constants import constants
from pynucastro.networks.rate_collection import RateCollection
from pynucastro.rates import ApproximateRate
from pynucastro.screening import get_screening_map


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
                c = r.reactants.count(ydot_i_nucleus)

                jac_str = r.jacobian_string_py(y_j_nucleus)
                if jac_str == "":
                    continue

                if c == 1:
                    rate_terms_str += f"{indent}   -{jac_str}\n"
                else:
                    rate_terms_str += f"{indent}   -{c}*{jac_str}\n"
            for r in self.nuclei_produced[ydot_i_nucleus]:
                c = r.products.count(ydot_i_nucleus)

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
        every set of reactants in our network, and updating the reaction
        rate values stored in the network.

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
        ostr += f"{indent}plasma_state = PlasmaState(T, rho, Y, Z)\n"

        if not self.do_screening:
            screening_map = []
        else:
            screening_map = get_screening_map(self.get_rates(),
                                              symmetric_screening=self.symmetric_screening)

        for i, scr in enumerate(screening_map):
            if not (scr.n1.dummy or scr.n2.dummy):
                # calculate the screening factor
                ostr += f"\n{indent}scn_fac = ScreenFactors({scr.n1.Z}, {scr.n1.A}, {scr.n2.Z}, {scr.n2.A})\n"
                ostr += f"{indent}scor = screen_func(plasma_state, scn_fac)\n"

            if scr.name == "He4_He4_He4":
                # we don't need to do anything here, but we want to avoid immediately applying the screening
                pass

            elif scr.name == "He4_He4_He4_dummy":
                # make sure the previous iteration was the first part of 3-alpha
                assert screening_map[i - 1].name == "He4_He4_He4"
                # handle the second part of the screening for 3-alpha
                ostr += f"{indent}scn_fac2 = ScreenFactors({scr.n1.Z}, {scr.n1.A}, {scr.n2.Z}, {scr.n2.A})\n"
                ostr += f"{indent}scor2 = screen_func(plasma_state, scn_fac2)\n"

                # there might be both the forward and reverse 3-alpha
                # if we are doing symmetric screening

                for r in scr.rates:
                    # use scor from the previous loop iteration
                    ostr += f"{indent}rate_eval.{r.fname} *= scor * scor2\n"
            else:
                # there might be several rates that have the same
                # reactants and therefore the same screening applies
                # -- handle them all now

                for r in scr.rates:
                    ostr += f"{indent}rate_eval.{r.fname} *= scor\n"

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
            return f"{indent}{r.fname}({', '.join(args)})\n"

        ostr = ""
        ostr += f"{indent}# reaclib rates\n"
        for r in self.reaclib_rates:
            ostr += format_rate_call(r)

        if self.derived_rates:
            ostr += f"\n{indent}# derived rates\n"
        for r in self.derived_rates:
            ostr += format_rate_call(r)

        if self.tabular_rates:
            ostr += f"\n{indent}# tabular rates\n"
        for r in self.tabular_rates:
            ostr += format_rate_call(r, use_tf=False)

        if self.custom_rates:
            ostr += f"\n{indent}# custom rates\n"
        for r in self.custom_rates:
            ostr += format_rate_call(r)

        ostr += "\n"

        # apply screening factors, if we're given a screening function
        ostr += f"{indent}if screen_func is not None:\n"
        ostr += self.screening_string(indent=indent + 4*" ")

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
        else:
            outfile = Path(outfile)
            of = outfile.open("w")

        indent = 4*" "

        of.write("import numba\n")
        of.write("import numpy as np\n")
        of.write("from scipy import constants\n")
        of.write("from numba.experimental import jitclass\n\n")
        of.write("from pynucastro.rates import TableIndex, TableInterpolator, TabularRate, Tfactors\n")
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
        of.write(f'{indent}'"enuc *= -1*constants.Avogadro\n")
        of.write(f'{indent}'"return enuc\n\n")

        # partition function data (if needed)

        nuclei_pfs = self.get_nuclei_needing_partition_functions()

        for n in nuclei_pfs:
            of.write(f"{n}_temp_array = np.array({list(n.partition_function.temperature/1.0e9)})\n")
            of.write(f"{n}_pf_array = np.array({list(n.partition_function.partition_function)})\n")
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
            of.write(f"{r.fname}_rate = TabularRate(rfile='{r.rfile}')\n")
            of.write(f"{r.fname}_info = ({r.fname}_rate.table_rhoy_lines,\n")
            of.write(f"                  {r.fname}_rate.table_temp_lines,\n")
            of.write(f"                  {r.fname}_rate.tabular_data_table)\n\n")

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

        if outfile is not None:
            of.close()

        # Copy any tables in the network to the current directory
        # if the table file cannot be found, print a warning and continue.
        try:
            odir = outfile.parent
        except AttributeError:
            odir = None

        for tr in self.tabular_rates:
            tdir = tr.rfile_path.parent
            if tdir != Path.cwd():
                tdat_file = tdir/tr.table_file
                if tdat_file.is_file():
                    shutil.copy(tdat_file, odir or Path.cwd())
                else:
                    warnings.warn(UserWarning(f'Table data file {tr.table_file} not found.'))
                rtoki_file = tdir/tr.rfile
                if rtoki_file.is_file():
                    shutil.copy(rtoki_file, odir or Path.cwd())
                else:
                    warnings.warn(UserWarning(f'Table metadata file {tr.rfile} not found.'))
