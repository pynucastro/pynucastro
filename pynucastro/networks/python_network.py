"""Support modules to write a pure python reaction network ODE
source"""

import sys

from pynucastro.networks.rate_collection import RateCollection
from pynucastro.rates.rate import ApproximateRate


class PythonNetwork(RateCollection):
    """A pure python reaction network."""

    def full_ydot_string(self, nucleus, indent=""):
        """construct the python form of dY(nucleus)/dt"""

        ostr = ""
        if not self.nuclei_consumed[nucleus] + self.nuclei_produced[nucleus]:
            # this captures an inert nucleus
            ostr += f"{indent}dYdt[j{nucleus}] = 0.0\n\n"
        else:
            ostr += f"{indent}dYdt[j{nucleus}] = (\n"
            for r in self.nuclei_consumed[nucleus]:
                c = r.reactants.count(nucleus)
                if c == 1:
                    ostr += f"{indent}   -{r.ydot_string_py()}\n"
                else:
                    ostr += f"{indent}   -{c}*{r.ydot_string_py()}\n"
            for r in self.nuclei_produced[nucleus]:
                c = r.products.count(nucleus)
                if c == 1:
                    ostr += f"{indent}   +{r.ydot_string_py()}\n"
                else:
                    ostr += f"{indent}   +{c}*{r.ydot_string_py()}\n"
            ostr += f"{indent}   )\n\n"

        return ostr

    def full_jacobian_element_string(self, ydot_i_nucleus, y_j_nucleus, indent=""):
        """return the Jacobian element dYdot(ydot_i_nucleus)/dY(y_j_nucleus)"""

        # this is the jac(i,j) string
        idx_str = f"jac[j{ydot_i_nucleus}, j{y_j_nucleus}]"

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
            else:
                ostr += rate_terms_str

            ostr += f"{indent}   )\n\n"

        return ostr

    def screening_string(self, indent=""):
        ostr = ""
        ostr += f"{indent}plasma_state = PlasmaState(T, rho, Y, Z)\n"

        screening_map = self.get_screening_map()
        for i, scr in enumerate(screening_map):
            if not (scr.n1.dummy or scr.n2.dummy):
                # calculate the screening factor
                ostr += f"\n{indent}scn_fac = ScreenFactors({scr.n1.Z}, {scr.n1.A}, {scr.n2.Z}, {scr.n2.A})\n"
                ostr += f"{indent}scor = screen_func(plasma_state, scn_fac)\n"

            if scr.name == "he4_he4_he4":
                # we don't need to do anything here, but we want to avoid immediately applying the screening
                pass

            elif scr.name == "he4_he4_he4_dummy":
                # make sure the previous iteration was the first part of 3-alpha
                assert screening_map[i - 1].name == "he4_he4_he4"
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
        """section for evaluating the rates and storing them in rate_eval"""
        ostr = ""
        ostr += f"{indent}# reaclib rates\n"
        for r in self.reaclib_rates:
            ostr += f"{indent}{r.fname}(rate_eval, tf)\n"

        if self.derived_rates:
            ostr += f"\n{indent}# derived rates\n"
        for r in self.derived_rates:
            ostr += f"{indent}{r.fname}(rate_eval, tf)\n"

        if self.tabular_rates:
            ostr += f"\n{indent}# tabular rates\n"
        for r in self.tabular_rates:
            ostr += f"{indent}{r.fname}(rate_eval, T, rho*ye(Y))\n"

        ostr += "\n"

        # apply screening factors, if we're given a screening function
        ostr += f"{indent}if screen_func is not None:\n"
        ostr += self.screening_string(indent=indent + 4*" ")

        if self.approx_rates:
            ostr += f"\n{indent}# approximate rates\n"
        for r in self.approx_rates:
            ostr += f"{indent}{r.fname}(rate_eval, tf)\n"

        return ostr

    def _write_network(self, outfile=None):
        """
        This is the actual RHS for the system of ODEs that
        this network describes.
        """
        # pylint: disable=arguments-differ
        if outfile is None:
            of = sys.stdout
        else:
            of = open(outfile, "w")

        indent = 4*" "

        of.write("import numba\n")
        of.write("import numpy as np\n")
        of.write("from numba.experimental import jitclass\n\n")
        of.write("from pynucastro.rates import Tfactors, _find_rate_file\n")
        of.write("from pynucastro.screening import PlasmaState, ScreenFactors\n\n")

        # integer keys

        for i, n in enumerate(self.unique_nuclei):
            of.write(f"j{n} = {i}\n")

        of.write(f"nnuc = {len(self.unique_nuclei)}\n\n")

        # nuclei properties

        of.write("A = np.zeros((nnuc), dtype=np.int32)\n\n")
        for n in self.unique_nuclei:
            of.write(f"A[j{n}] = {n.A}\n")

        of.write("\n")

        of.write("Z = np.zeros((nnuc), dtype=np.int32)\n\n")
        for n in self.unique_nuclei:
            of.write(f"Z[j{n}] = {n.Z}\n")

        of.write("\n")

        of.write("names = []\n")
        for n in self.unique_nuclei:
            of.write(f"names.append(\"{n.short_spec_name}\")\n")

        of.write("\n")

        # partition function data (if needed)

        nuclei_pfs = self.get_nuclei_needing_partition_functions()

        if nuclei_pfs:
            for n in nuclei_pfs:
                if n.partition_function:
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
        for r in self.tabular_rates:

            of.write(f"# load data for {r.rid}\n")
            of.write(f"{r.fname}_table_path = _find_rate_file('{r.table_file}')\n")
            of.write("t_data2d = []\n")
            of.write(f"with open({r.fname}_table_path) as tabular_file:\n")
            of.write(f'{indent}'"for i, line in enumerate(tabular_file):\n")
            of.write(f'{indent*2}'f"if i < {r.table_header_lines}:\n")
            of.write(f'{indent*3}'"continue\n")
            of.write(f'{indent*2}'"line = line.strip()\n")
            of.write(f'{indent*2}'"if not line:\n")
            of.write(f'{indent*3}'"continue\n")
            of.write(f'{indent*2}'"t_data2d.append(line.split())\n")
            of.write(f"{r.fname}_data = np.array(t_data2d, dtype=float)\n\n")

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
