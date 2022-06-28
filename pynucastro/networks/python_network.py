"""Support modules to write a pure python reaction network ODE
source"""

import sys

from pynucastro.networks import RateCollection
from pynucastro.rates.rate import ApproximateRate


class PythonNetwork(RateCollection):
    """A pure python reaction network."""

    def approx_function_string(self, rate):
        """
        Return a string containing python function that computes the
        approximate rate
        """

        if not rate.approx_type == "ap_pg":
            raise NotImplementedError("don't know how to work with this approximation")

        string = ""
        string += "@numba.njit()\n"
        string += f"def {rate.fname}(tf):\n"

        if not rate.is_reverse:

            # first we need to get all of the rates that make this up
            string += f"    r_ag = {rate.primary_rate.fname}(tf)\n"
            string += f"    r_ap = {rate.secondary_rates[0].fname}(tf)\n"
            string += f"    r_pg = {rate.secondary_rates[1].fname}(tf)\n"
            string += f"    r_pa = {rate.secondary_reverse[1].fname}(tf)\n"

            # now the approximation
            string += "    rate = r_ag + r_ap * r_pg / (r_pg + r_pa)\n"

        else:

            # first we need to get all of the rates that make this up
            string += f"    r_ga = {rate.primary_reverse.fname}(tf)\n"
            string += f"    r_pa = {rate.secondary_reverse[1].fname}(tf)\n"
            string += f"    r_gp = {rate.secondary_reverse[0].fname}(tf)\n"
            string += f"    r_pg = {rate.secondary_rates[1].fname}(tf)\n"

            # now the approximation
            string += "    rate = r_ga + r_pa * r_gp / (r_pg + r_pa)\n"

        string += "    return rate\n\n"
        return string

    def ydot_string(self, rate):
        """
        Return a string containing the term in a dY/dt equation
        in a reaction network corresponding to this rate.
        """

        # composition dependence
        Y_string = ""
        for n, r in enumerate(sorted(set(rate.reactants))):
            c = rate.reactants.count(r)
            if c > 1:
                Y_string += f"Y[j{r}]**{c}"
            else:
                Y_string += f"Y[j{r}]"

            if n < len(set(rate.reactants))-1:
                Y_string += "*"

        # density dependence
        if rate.dens_exp == 0:
            dens_string = ""
        elif rate.dens_exp == 1:
            dens_string = "rho*"
        else:
            dens_string = f"rho**{rate.dens_exp}*"

        # electron fraction dependence
        if (rate.weak_type == 'electron_capture' and not rate.tabular):
            y_e_string = 'ye(Y)*'
        else:
            y_e_string = ''

        # prefactor
        if not rate.prefactor == 1.0:
            prefactor_string = f"{rate.prefactor:1.14e}*"
        else:
            prefactor_string = ""

        return "{}{}{}{}*lambda_{}".format(prefactor_string, dens_string,
                                           y_e_string, Y_string, rate.fname)

    def full_ydot_string(self, nucleus, indent=""):

        ostr = ""
        if not self.nuclei_consumed[nucleus] + self.nuclei_produced[nucleus]:
            ostr += f"{indent}dYdt[j{nucleus}] = 0.0\n\n"
        else:
            ostr += f"{indent}dYdt[j{nucleus}] = (\n"
            for r in self.nuclei_consumed[nucleus]:
                c = r.reactants.count(nucleus)
                if c == 1:
                    ostr += f"{indent}   -{self.ydot_string(r)}\n"
                else:
                    ostr += f"{indent}   -{c}*{self.ydot_string(r)}\n"
            for r in self.nuclei_produced[nucleus]:
                c = r.products.count(nucleus)
                if c == 1:
                    ostr += f"{indent}   +{self.ydot_string(r)}\n"
                else:
                    ostr += f"{indent}   +{c}*{self.ydot_string(r)}\n"
            ostr += f"{indent}   )\n\n"

        return ostr

    def jacobian_string(self, rate, ydot_j, y_i):
        """
        Return a string containing the term in a jacobian matrix
        in a reaction network corresponding to this rate.

        Returns the derivative of the j-th YDOT wrt. the i-th Y
        If the derivative is zero, returns the empty string ''

        ydot_j and y_i are objects of the class ``Nucleus``.
        """
        if (ydot_j not in rate.reactants and ydot_j not in rate.products) or \
           y_i not in rate.reactants:
            return ''

        # composition dependence
        Y_string = ""
        for n, r in enumerate(sorted(set(rate.reactants))):
            c = rate.reactants.count(r)
            if y_i == r:
                if c == 1:
                    continue
                if 0 < n < len(set(rate.reactants))-1:
                    Y_string += "*"
                if c > 2:
                    Y_string += f"{c}*Y[j{r}]**{c-1}"
                elif c == 2:
                    Y_string += f"2*Y[j{r}]"
            else:
                if 0 < n < len(set(rate.reactants))-1:
                    Y_string += "*"
                if c > 1:
                    Y_string += f"Y[j{r}]**{c}"
                else:
                    Y_string += f"Y[j{r}]"

        # density dependence
        if rate.dens_exp == 0:
            dens_string = ""
        elif rate.dens_exp == 1:
            dens_string = "rho*"
        else:
            dens_string = f"rho**{rate.dens_exp}*"

        # electron fraction dependence
        if (rate.weak_type == 'electron_capture' and not rate.tabular):
            y_e_string = 'ye(Y)*'
        else:
            y_e_string = ''

        # prefactor
        if not rate.prefactor == 1.0:
            prefactor_string = f"{rate.prefactor:1.14e}*"
        else:
            prefactor_string = ""

        if Y_string == "" and dens_string == "" and prefactor_string == "":
            rstring = "{}{}{}lambda_{}"
        else:
            rstring = "{}{}{}{}*lambda_{}"
        return rstring.format(prefactor_string, dens_string,
                              y_e_string, Y_string, rate.fname)

    def _write_network(self, outfile=None):
        """
        This is the actual RHS for the system of ODEs that
        this network describes.
        """
        if outfile is None:
            of = sys.stdout
        else:
            try:
                of = open(outfile, "w")
            except IOError:
                print(f"unable to open {outfile}")
                raise

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
                of.write(self.approx_function_string(r))
            else:
                if r in _rate_func_written:
                    continue
                of.write(r.function_string_py())
                _rate_func_written.append(r)

        of.write("def rhs(t, Y, rho, T):\n")
        of.write(f"{indent}return rhs_eq(t, Y, rho, T)\n\n")

        of.write("@numba.njit()\n")
        of.write("def rhs_eq(t, Y, rho, T):\n\n")
        # integer keys
        for i, n in enumerate(self.unique_nuclei):
            of.write(f"{indent}j{n} = {i}\n")

        of.write(f"{indent}nnuc = {len(self.unique_nuclei)}\n\n")

        # get the rates
        of.write(f"{indent}tf = Tfactors(T)\n\n")
        for r in self.rates:
            of.write(f"{indent}lambda_{r.fname} = {r.fname}(tf)\n")

        of.write("\n")

        of.write(f"{indent}dYdt = np.zeros((nnuc), dtype=np.float64)\n\n")

        # now make the RHSs
        for n in self.unique_nuclei:
            of.write(self.full_ydot_string(n, indent=indent))

        of.write(f"{indent}return dYdt\n")
