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

        of.write("import numpy as np\n")
        of.write("from pynucastro.rates import Tfactors\n")
        of.write("import numba\n\n")

        # integer keys
        for i, n in enumerate(self.unique_nuclei):
            of.write(f"j{n} = {i}\n")

        of.write(f"nnuc = {len(self.unique_nuclei)}\n\n")

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

        indent = 4*" "

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

        of.write("def rhs(t, Y, rho, T):\n")
        of.write(f"{indent}return rhs_eq(t, Y, rho, T)\n\n")

        of.write("@numba.njit()\n")
        of.write("def rhs_eq(t, Y, rho, T):\n\n")

        # get the rates
        of.write(f"{indent}tf = Tfactors(T)\n\n")
        for r in self.rates:
            of.write(f"{indent}lambda_{r.fname} = {r.fname}(tf)\n")

        of.write("\n")

        of.write(f"{indent}dYdt = np.zeros((nnuc), dtype=np.float64)\n\n")

        # now make the RHSs
        for n in self.unique_nuclei:
            of.write(self.full_ydot_string(n, indent=indent))

        of.write(f"{indent}return dYdt\n\n")

        # the jacobian() function

        of.write("def jacobian(t, Y, rho, T):\n")
        of.write(f"{indent}return jacobian_eq(t, Y, rho, T)\n\n")

        of.write("@numba.njit()\n")
        of.write("def jacobian_eq(t, Y, rho, T):\n\n")

        # get the rates
        of.write(f"{indent}tf = Tfactors(T)\n\n")
        for r in self.rates:
            of.write(f"{indent}lambda_{r.fname} = {r.fname}(tf)\n")

        of.write("\n")

        of.write(f"{indent}jac = np.zeros((nnuc, nnuc), dtype=np.float64)\n\n")

        # now fill each Jacobian element
        for n_i in self.unique_nuclei:
            for n_j in self.unique_nuclei:
                of.write(self.full_jacobian_element_string(n_i, n_j, indent=indent))

        of.write(f"{indent}return jac\n")