"""Support modules to write a pure python reaction network ODE
source"""

import sys

from pynucastro.networks import RateCollection
from pynucastro.rates.rate import ApproximateRate


class PythonNetwork(RateCollection):
    """A pure python reaction network."""

    def rate_string(self, rate, indent=0, prefix="rate"):
        """
        Return the functional form of rate as a function of
        the temperature (as Tfactors)

        rate is an object of class Rate
        """

        tstring = f"# {rate.rid}\n"
        tstring += f"{prefix} = 0.0\n\n"

        for s in rate.sets:
            tstring += f"# {s.labelprops[0:5]}\n"
            tstring += f"{s.set_string(prefix=prefix, plus_equal=True)}\n"

        string = ""
        for t in tstring.split("\n"):
            string += indent*" " + t + "\n"
        return string

    def function_string(self, rate):
        """
        Return a string containing python function that computes the
        rate
        """

        string = ""
        string += "@numba.njit()\n"
        string += f"def {rate.fname}(tf):\n"
        string += f"{self.rate_string(rate, indent=4)}"
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
                    of.write(self.function_string(cr))
                    _rate_func_written.append(cr)
            else:
                if r in _rate_func_written:
                    continue
                of.write(self.function_string(r))
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
            of.write(f"{indent}dYdt[j{n}] = (\n")
            for r in self.nuclei_consumed[n]:
                c = r.reactants.count(n)
                if c == 1:
                    of.write(f"{indent}   -{self.ydot_string(r)}\n")
                else:
                    of.write(f"{indent}   -{c}*{self.ydot_string(r)}\n")
            for r in self.nuclei_produced[n]:
                c = r.products.count(n)
                if c == 1:
                    of.write(f"{indent}   +{self.ydot_string(r)}\n")
                else:
                    of.write(f"{indent}   +{c}*{self.ydot_string(r)}\n")
            of.write(f"{indent}   )\n\n")

        of.write(f"{indent}return dYdt\n")
