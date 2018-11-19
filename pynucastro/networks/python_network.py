"""Support modules to write a pure python reaction network ODE
source"""

from __future__ import print_function

import sys

from pynucastro.networks import RateCollection


class PythonNetwork(RateCollection):
    """A pure python reaction network."""

    def rate_string(self, rate, indent=0, prefix="rate"):
        """
        Return the functional form of rate as a function of
        the temperature (as Tfactors)

        rate is an object of class Rate
        """

        tstring = "# {}\n".format(rate.string)
        tstring += "{} = 0.0\n\n".format(prefix)

        for s in rate.sets:
            tstring += "# {}\n".format(s.labelprops[0:5])
            tstring += "{}\n".format(s.set_string(prefix=prefix, plus_equal=True))

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
        string += "def {}(tf):\n".format(rate.fname)
        string += "{}".format(self.rate_string(rate, indent=4))
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
                Y_string += "Y[i{}]**{}".format(r, c)
            else:
                Y_string += "Y[i{}]".format(r)

            if n < len(set(rate.reactants))-1:
                Y_string += "*"

        # density dependence
        if rate.dens_exp == 0:
            dens_string = ""
        elif rate.dens_exp == 1:
            dens_string = "rho*"
        else:
            dens_string = "rho**{}*".format(rate.dens_exp)

        # electron fraction dependence
        if (rate.weak_type == 'electron_capture' and not rate.tabular):
            y_e_string = 'ye(Y)*'
        else:
            y_e_string = ''

        # prefactor
        if not rate.prefactor == 1.0:
            prefactor_string = "{:1.14e}*".format(rate.prefactor)
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
                if n > 0 and n < len(set(rate.reactants))-1:
                    Y_string += "*"
                if c > 2:
                    Y_string += "{}*Y[i{}]**{}".format(c, r, c-1)
                elif c == 2:
                    Y_string += "2*Y[i{}]".format(r)
            else:
                if n > 0 and n < len(set(rate.reactants))-1:
                    Y_string += "*"
                if c > 1:
                    Y_string += "Y[i{}]**{}".format(r, c)
                else:
                    Y_string += "Y[i{}]".format(r)

        # density dependence
        if rate.dens_exp == 0:
            dens_string = ""
        elif rate.dens_exp == 1:
            dens_string = "rho*"
        else:
            dens_string = "rho**{}*".format(rate.dens_exp)

        # electron fraction dependence
        if (rate.weak_type == 'electron_capture' and not rate.tabular):
            y_e_string = 'ye(Y)*'
        else:
            y_e_string = ''

        # prefactor
        if not rate.prefactor == 1.0:
            prefactor_string = "{:1.14e}*".format(rate.prefactor)
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
            except:
                raise

        of.write("import numpy as np\n")
        of.write("from pynucastro.rates import Tfactors\n")
        of.write("import numba\n\n")

        # integer keys
        for i, n in enumerate(self.unique_nuclei):
            of.write("i{} = {}\n".format(n, i))

        of.write("nnuc = {}\n\n".format(len(self.unique_nuclei)))

        of.write("A = np.zeros((nnuc), dtype=np.int32)\n\n")
        for n in self.unique_nuclei:
            of.write("A[i{}] = {}\n".format(n, n.A))

        of.write("\n")

        of.write("Z = np.zeros((nnuc), dtype=np.int32)\n\n")
        for n in self.unique_nuclei:
            of.write("Z[i{}] = {}\n".format(n, n.Z))

        of.write("\n")

        indent = 4*" "

        of.write("@numba.njit()\n")

        of.write("def ye(Y):\n")
        of.write("{}return np.sum(Z * Y)/np.sum(A * Y)\n\n".format(indent))

        for r in self.rates:
            of.write(self.function_string(r))

        of.write("def rhs(t, Y, rho, T):\n")
        of.write("{}return rhs_eq(t, Y, rho, T)\n\n".format(indent))

        of.write("@numba.njit()\n")
        of.write("def rhs_eq(t, Y, rho, T):\n\n")
        # integer keys
        for i, n in enumerate(self.unique_nuclei):
            of.write("{}i{} = {}\n".format(indent, n, i))

        of.write("{}nnuc = {}\n\n".format(indent, len(self.unique_nuclei)))


        # get the rates
        of.write("{}tf = Tfactors(T)\n\n".format(indent))
        for r in self.rates:
            of.write("{}lambda_{} = {}(tf)\n".format(indent, r.fname, r.fname))

        of.write("\n")

        of.write("{}dYdt = np.zeros((nnuc), dtype=np.float64)\n\n".format(indent))

        # now make the RHSs
        for n in self.unique_nuclei:
            of.write("{}dYdt[i{}] = (\n".format(indent, n))
            for r in self.nuclei_consumed[n]:
                c = r.reactants.count(n)
                if c == 1:
                    of.write("{}   -{}\n".format(indent, self.ydot_string(r)))
                else:
                    of.write("{}   -{}*{}\n".format(indent, c, self.ydot_string(r)))
            for r in self.nuclei_produced[n]:
                c = r.products.count(n)
                if c == 1:
                    of.write("{}   +{}\n".format(indent, self.ydot_string(r)))
                else:
                    of.write("{}   +{}*{}\n".format(indent, c, self.ydot_string(r)))
            of.write("{}   )\n\n".format(indent))

        of.write("{}return dYdt\n".format(indent))
