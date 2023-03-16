"""
Classes and methods to interface with files storing rate data.
"""
import io
import os
from collections import Counter

import matplotlib.pyplot as plt
import numpy as np
from scipy.constants import physical_constants

try:
    import numba
    from numba.experimental import jitclass
except ImportError:
    numba = None
    import functools

    # no-op jitclass placeholder
    def jitclass(cls_or_spec=None, spec=None):
        if (cls_or_spec is not None and
            spec is None and
                not isinstance(cls_or_spec, type)):
            # Used like
            # @jitclass([("x", intp)])
            # class Foo:
            #     ...
            spec = cls_or_spec
            cls_or_spec = None

        def wrap(cls):
            # this copies the function name and docstring to the wrapper function
            @functools.wraps(cls)
            def wrapper(*args, **kwargs):
                return cls(*args, **kwargs)
            return wrapper

        if cls_or_spec is None:
            return wrap
        return wrap(cls_or_spec)


from pynucastro.nucdata import Nucleus

amu_mev, _, _ = physical_constants['atomic mass constant energy equivalent in MeV']
hbar, _, _ = physical_constants['reduced Planck constant']
amu, _, _ = physical_constants['atomic mass constant']
k_B_mev_k, _, _ = physical_constants['Boltzmann constant in eV/K']
k_B_mev_k /= 1.0e6
k_B, _, _ = physical_constants['Boltzmann constant']
N_a, _, _ = physical_constants['Avogadro constant']

_pynucastro_dir = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
_pynucastro_rates_dir = os.path.join(_pynucastro_dir, 'library')
_pynucastro_tabular_dir = os.path.join(_pynucastro_rates_dir, 'tabular')


def load_rate(rfile=None, rfile_path=None):

    try:
        rate = TabularRate(rfile=rfile, rfile_path=rfile_path)
    except AssertionError:
        try:
            rate = ReacLibRate(rfile=rfile, rfile_path=rfile_path)
        except AssertionError:
            raise

    return rate


def _find_rate_file(ratename):
    """locate the Reaclib or tabular rate or library file given its name.  Return
    None if the file cannot be located, otherwise return its path."""

    # check to see if the rate file is in the working dir or
    # is already the full path
    x = ratename
    if os.path.isfile(x):
        return os.path.realpath(x)

    # check to see if the rate file is in pynucastro/library
    x = os.path.join(_pynucastro_rates_dir, ratename)
    if os.path.isfile(x):
        return os.path.realpath(x)

    # check to see if the rate file is in pynucastro/library/tabular
    x = os.path.join(_pynucastro_tabular_dir, ratename)
    if os.path.isfile(x):
        return os.path.realpath(x)

    # notify user we can't find the file
    raise Exception(f'File {ratename!r} not found in the working directory, {_pynucastro_rates_dir}, or {_pynucastro_tabular_dir}')


if numba is not None:
    Tfactor_spec = [
        ('T9', numba.float64),
        ('T9i', numba.float64),
        ('T913', numba.float64),
        ('T913i', numba.float64),
        ('T953', numba.float64),
        ('lnT9', numba.float64)
    ]
else:
    Tfactor_spec = []


@jitclass(Tfactor_spec)
class Tfactors:
    """ precompute temperature factors for speed

    :param float T: input temperature (Kelvin)
    :var T9:    T / 1.e9 K
    :var T9i:   1.0 / T9
    :var T913i  1.0 / T9 ** (1/3)
    :var T913   T9 ** (1/3)
    :var T953   T9 ** (5/3)
    :var lnT9   log(T9)
    """

    def __init__(self, T):
        """ return the Tfactors object.  Here, T is temperature in Kelvin """
        self.T9 = T/1.e9
        self.T9i = 1.0/self.T9
        self.T913i = self.T9i**(1./3.)
        self.T913 = self.T9**(1./3.)
        self.T953 = self.T9**(5./3.)
        self.lnT9 = np.log(self.T9)


class SingleSet:
    """ a set in Reaclib is one piece of a rate, in the form

        lambda = exp[ a_0 + sum_{i=1}^5  a_i T_9**(2i-5)/3  + a_6 log T_9]

    A single rate in Reaclib can be composed of multiple sets

    :param a: the coefficients of the exponential fit
    :param labelprops: a collection of flags that classify a ReacLib rate

    """

    def __init__(self, a, labelprops):
        """here a is iterable (e.g., list or numpy array), storing the
           coefficients, a0, ..., a6

        """
        self.a = a
        self.labelprops = labelprops
        self.label = None
        self.resonant = None
        self.weak = None
        self.reverse = None

        self._update_label_properties()

    def _update_label_properties(self):
        """ Set label and flags indicating Set is resonant,
            weak, or reverse. """
        assert isinstance(self.labelprops, str)
        assert len(self.labelprops) == 6

        self.label = self.labelprops[0:4]
        self.resonant = self.labelprops[4] == 'r'
        self.weak = self.labelprops[4] == 'w'
        self.reverse = self.labelprops[5] == 'v'

    def __eq__(self, other):
        """ Determine whether two SingleSet objects are equal to each other. """
        x = True

        for ai, aj in zip(self.a, other.a):
            x = x and (ai == aj)

        x = x and (self.label == other.label)
        x = x and (self.resonant == other.resonant)
        x = x and (self.weak == other.weak)
        x = x and (self.reverse == other.reverse)
        return x

    def f(self):
        """ return a function for rate(tf) where tf is a Tfactors
        object """
        return lambda tf: np.exp(self.a[0] +
                                 self.a[1]*tf.T9i +
                                 self.a[2]*tf.T913i +
                                 self.a[3]*tf.T913 +
                                 self.a[4]*tf.T9 +
                                 self.a[5]*tf.T953 +
                                 self.a[6]*tf.lnT9)

    def dfdT(self):
        """ return a function for this dratedT(tf), where tf is a
        Tfactors object """

        # we have lambda = exp(f(T_9))
        # so dlambda/dT9 = lambda * df/dT9
        # and dlambda/dT = dlambda/dT9 / 1.e9

        return lambda tf: self.f()(tf) * (-self.a[1] * tf.T9i * tf.T9i +
                                          -(1./3.) * self.a[2] * tf.T913i * tf.T9i +
                                          (1./3.) * self.a[3] * tf.T913i * tf.T913i +
                                          self.a[4] +
                                          (5./3.) * self.a[5] * tf.T913 * tf.T913 +
                                          self.a[6] * tf.T9i) / 1.e9

    def set_string_py(self, prefix="set", plus_equal=False):
        """
        return a string containing the python code for this set
        """
        if plus_equal:
            string = f"{prefix} += np.exp( "
        else:
            string = f"{prefix} = np.exp( "
        string += f" {self.a[0]}"
        if not self.a[1] == 0.0:
            string += f" + {self.a[1]}*tf.T9i"
        if not self.a[2] == 0.0:
            string += f" + {self.a[2]}*tf.T913i"
        if not self.a[3] == 0.0:
            string += f" + {self.a[3]}*tf.T913"
        if not (self.a[4] == 0.0 and self.a[5] == 0.0 and self.a[6] == 0.0):
            indent = len(prefix)*" "
            string += f"\n{indent}         "
        if not self.a[4] == 0.0:
            string += f" + {self.a[4]}*tf.T9"
        if not self.a[5] == 0.0:
            string += f" + {self.a[5]}*tf.T953"
        if not self.a[6] == 0.0:
            string += f" + {self.a[6]}*tf.lnT9"
        string += ")"
        return string

    def set_string_cxx(self, prefix="set", plus_equal=False, with_exp=True):
        """
        return a string containing the C++ code for this set
        """
        if plus_equal:
            string = f"{prefix} += "
        else:
            string = f"{prefix} = "
        if with_exp:
            string += "std::exp( "
        string += f" {self.a[0]}"
        if not self.a[1] == 0.0:
            string += f" + {self.a[1]} * tfactors.T9i"
        if not self.a[2] == 0.0:
            string += f" + {self.a[2]} * tfactors.T913i"
        if not self.a[3] == 0.0:
            string += f" + {self.a[3]} * tfactors.T913"
        if not (self.a[4] == 0.0 and self.a[5] == 0.0 and self.a[6] == 0.0):
            indent = len(prefix)*" "
            string += f"\n{indent}         "
        if not self.a[4] == 0.0:
            string += f" + {self.a[4]} * tfactors.T9"
        if not self.a[5] == 0.0:
            string += f" + {self.a[5]} * tfactors.T953"
        if not self.a[6] == 0.0:
            string += f" + {self.a[6]} * tfactors.lnT9"
        if with_exp:
            string += ");"
        else:
            string += ";"
        return string

    def dln_set_string_dT9_cxx(self, prefix="dset_dT", plus_equal=False):
        """
        return a string containing the C++ code for d/dT9 ln(set)
        """
        if plus_equal:
            string = f"{prefix} += "
        else:
            string = f"{prefix} = "

        if all(q == 0.0 for q in self.a[1:]):
            string += "0.0;"
            return string

        if not self.a[1] == 0.0:
            string += f" {-self.a[1]} * tfactors.T9i * tfactors.T9i"
        if not self.a[2] == 0.0:
            string += f" + -(1.0/3.0) * {self.a[2]} * tfactors.T943i"
        if not self.a[3] == 0.0:
            string += f" + (1.0/3.0) * {self.a[3]} * tfactors.T923i"
        if not (self.a[4] == 0.0 and self.a[5] == 0.0 and self.a[6] == 0.0):
            indent = len(prefix)*" "
            string += f"\n{indent}         "
        if not self.a[4] == 0.0:
            string += f" + {self.a[4]}"
        if not self.a[5] == 0.0:
            string += f" + (5.0/3.0) * {self.a[5]} * tfactors.T923"
        if not self.a[6] == 0.0:
            string += f" + {self.a[6]} * tfactors.T9i"
        string += ";"
        return string


class Rate:
    """The base reaction rate class.  Most rate types will subclass
    this and extend to their particular format.

    """
    def __init__(self, reactants=None, products=None, Q=None, weak_type=""):
        """ rfile can be either a string specifying the path to a rate file or
        an io.StringIO object from which to read rate information. """

        self.fname = None

        if reactants:
            self.reactants = reactants
        else:
            self.reactants = []

        if products:
            self.products = products
        else:
            self.products = []

        self.label = "generic"

        self.Q = Q

        self.weak_type = weak_type

        self._set_rhs_properties()
        self._set_screening()
        self._set_print_representation()

        self.tabular = False

    def __repr__(self):
        return self.string

    def __hash__(self):
        return hash(self.__repr__())

    def __eq__(self, other):
        """ Determine whether two Rate objects are equal.
        They are equal if they contain identical reactants and products and
        if they contain the same SingleSet sets and if their chapters are equal."""
        x = True

        x = x and (self.reactants == other.reactants)
        x = x and (self.products == other.products)
        return x

    def __lt__(self, other):
        """sort such that lightest reactants come first, and then look at products"""

        # this sort will make two nuclei with the same A be in order of Z
        # (assuming there are no nuclei with A > 999
        # we want to compare based on the heaviest first, so we reverse

        self_react_sorted = sorted(self.reactants, key=lambda x: 1000*x.A + x.Z, reverse=True)
        other_react_sorted = sorted(other.reactants, key=lambda x: 1000*x.A + x.Z, reverse=True)

        if self_react_sorted != other_react_sorted:
            # reactants are different, so now we can check them
            for srn, orn in zip(self_react_sorted, other_react_sorted):
                if not srn == orn:
                    return srn < orn
        else:
            # reactants are the same, so consider products
            self_prod_sorted = sorted(self.products, key=lambda x: 1000*x.A + x.Z, reverse=True)
            other_prod_sorted = sorted(other.products, key=lambda x: 1000*x.A + x.Z, reverse=True)

            for spn, opn in zip(self_prod_sorted, other_prod_sorted):
                if not spn == opn:
                    return spn < opn

        # if we made it here, then the rates are the same
        return True

    def _set_q(self):
        """set the Q value of the reaction (in MeV)"""

        # from the binding energy of the nuclei, Q = -B_reactants + B_products
        # but note that nucbind is the binding energy *per* nucleon, so we need
        # to multiply by the number of nucleons

        self.Q = 0
        for n in self.reactants:
            self.Q += -n.A * n.nucbind
        for n in self.products:
            self.Q += n.A * n.nucbind

    def _set_print_representation(self):

        # string is output to the terminal, rid is used as a dict key,
        # and pretty_string is latex

        # some rates will have no nuclei particles (e.g. gamma) on the left or
        # right -- we'll try to infer those here
        lhs_other = []
        rhs_other = []

        self.string = ""
        self.rid = ""
        self.pretty_string = r"$"

        # put p, n, and alpha second
        treactants = []
        for n in self.reactants:
            if n.raw not in ["p", "he4", "n"]:
                treactants.insert(0, n)
            else:
                treactants.append(n)

        # figure out if there are any non-nuclei present
        # for the moment, we just handle strong rates

        # there should be the same number of protons on each side and
        # the same number of neutrons on each side

        strong_test = sum(n.Z for n in self.reactants) == sum(n.Z for n in self.products) and \
                      sum(n.A for n in self.reactants) == sum(n.A for n in self.products)

        if strong_test:
            if len(self.products) == 1:
                rhs_other.append("gamma")
        else:
            # this is a weak rate

            if self.weak_type == "electron_capture":

                # we assume that all the tabular rates are electron capture for now

                # we expect an electron on the left -- let's make sure
                # the charge on the left should be +1 the charge on the right
                assert sum(n.Z for n in self.reactants) == sum(n.Z for n in self.products) + 1

                lhs_other.append("e-")
                rhs_other.append("nu")

            elif self.weak_type == "beta_decay":
                # we expect an electron on the right
                assert sum(n.Z for n in self.reactants) + 1 == sum(n.Z for n in self.products)

                rhs_other.append("e-")
                rhs_other.append("nubar")

            elif "_pos_" in self.weak_type:

                # we expect a positron on the right -- let's make sure
                assert sum(n.Z for n in self.reactants) == sum(n.Z for n in self.products) + 1

                rhs_other.append("e+")
                rhs_other.append("nu")

            elif "_neg_" in self.weak_type:

                # we expect an electron on the right -- let's make sure
                assert sum(n.Z for n in self.reactants) + 1 == sum(n.Z for n in self.products)

                rhs_other.append("e-")
                rhs_other.append("nubar")

            else:

                # we need to figure out what the rate is.  We'll assume that it is
                # not an electron capture

                if sum(n.Z for n in self.reactants) == sum(n.Z for n in self.products) + 1:
                    rhs_other.append("e+")
                    rhs_other.append("nu")

                elif sum(n.Z for n in self.reactants) + 1 == sum(n.Z for n in self.products):

                    rhs_other.append("e-")
                    rhs_other.append("nubar")

        for n, r in enumerate(treactants):
            self.string += f"{r.c()}"
            self.rid += f"{r}"
            self.pretty_string += fr"{r.pretty}"
            if not n == len(self.reactants)-1:
                self.string += " + "
                self.rid += " + "
                self.pretty_string += r" + "

        if lhs_other:
            for o in lhs_other:
                if o == "e-":
                    self.string += " + e⁻"
                    self.pretty_string += r" + \mathrm{e}^-"

        self.string += " ⟶ "
        self.rid += " --> "
        self.pretty_string += r" \rightarrow "

        for n, p in enumerate(self.products):
            self.string += f"{p.c()}"
            self.rid += f"{p}"
            self.pretty_string += fr"{p.pretty}"
            if not n == len(self.products)-1:
                self.string += " + "
                self.rid += " + "
                self.pretty_string += r" + "

        if rhs_other:
            for o in rhs_other:
                if o == "gamma":
                    self.string += " + 𝛾"
                    self.pretty_string += r"+ \gamma"
                elif o == "nu":
                    self.string += " + 𝜈"
                    self.pretty_string += r"+ \nu_e"
                elif o == "nubar":
                    self.string += " + 𝜈"
                    self.pretty_string += r"+ \bar{\nu}_e"
                if o == "e-":
                    self.string += " + e⁻"
                    self.pretty_string += r" + \mathrm{e}^-"
                if o == "e+":
                    self.string += " + e⁺"
                    self.pretty_string += r" + \mathrm{e}^+"

        self.pretty_string += r"$"

    def _set_rhs_properties(self):
        """ compute statistical prefactor and density exponent from the reactants. """
        self.prefactor = 1.0  # this is 1/2 for rates like a + a (double counting)
        self.inv_prefactor = 1
        for r in set(self.reactants):
            self.inv_prefactor = self.inv_prefactor * np.math.factorial(self.reactants.count(r))
        self.prefactor = self.prefactor/float(self.inv_prefactor)
        self.dens_exp = len(self.reactants)-1
        if self.weak_type == 'electron_capture':
            self.dens_exp = self.dens_exp + 1

    def _set_screening(self):
        """ determine if this rate is eligible for screening and the nuclei to use. """
        # Tells if this rate is eligible for screening
        # using screenz.f90 provided by StarKiller Microphysics.
        # If not eligible for screening, set to None
        # If eligible for screening, then
        # Rate.ion_screen is a 2-element (3 for 3-alpha) list of Nucleus objects for screening
        self.ion_screen = []
        nucz = [q for q in self.reactants if q.Z != 0]
        if len(nucz) > 1:
            nucz.sort(key=lambda x: x.Z)
            self.ion_screen = []
            self.ion_screen.append(nucz[0])
            self.ion_screen.append(nucz[1])
            if len(nucz) == 3:
                self.ion_screen.append(nucz[2])

        # if the rate is a reverse rate (defined as Q < 0), then we
        # might actually want to compute the screening based on the
        # reactants of the forward rate that was used in the detailed
        # balance.  Rate.symmetric_screen is what should be used in
        # the screening in this case
        self.symmetric_screen = []
        if self.Q < 0:
            nucz = [q for q in self.products if q.Z != 0]
            if len(nucz) > 1:
                nucz.sort(key=lambda x: x.Z)
                self.symmetric_screen = []
                self.symmetric_screen.append(nucz[0])
                self.symmetric_screen.append(nucz[1])
                if len(nucz) == 3:
                    self.symmetric_screen.append(nucz[2])
        else:
            self.symmetric_screen = self.ion_screen

    def get_rate_id(self):
        """ Get an identifying string for this rate."""
        return f'{self.rid} <{self.label.strip()}>'

    def heaviest(self):
        """
        Return the heaviest nuclide in this Rate.

        If two nuclei are tied in mass number, return the one with the
        lowest atomic number.
        """
        nuc = self.reactants[0]
        for n in self.reactants + self.products:
            if n.A > nuc.A or (n.A == nuc.A and n.Z < nuc.Z):
                nuc = n
        return nuc

    def lightest(self):
        """
        Return the lightest nuclide in this Rate.

        If two nuclei are tied in mass number, return the one with the
        highest atomic number.
        """
        nuc = self.reactants[0]
        for n in self.reactants + self.products:
            if n.A < nuc.A or (n.A == nuc.A and n.Z > nuc.Z):
                nuc = n
        return nuc

    def ydot_string_py(self):
        """
        Return a string containing the term in a dY/dt equation
        in a reaction network corresponding to this rate.
        """

        # composition dependence
        Y_string = ""
        for n, r in enumerate(sorted(set(self.reactants))):
            c = self.reactants.count(r)
            if c > 1:
                Y_string += f"Y[j{r}]**{c}"
            else:
                Y_string += f"Y[j{r}]"

            if n < len(set(self.reactants))-1:
                Y_string += "*"

        # density dependence
        if self.dens_exp == 0:
            dens_string = ""
        elif self.dens_exp == 1:
            dens_string = "rho*"
        else:
            dens_string = f"rho**{self.dens_exp}*"

        # electron fraction dependence
        if self.weak_type == 'electron_capture' and not self.tabular:
            y_e_string = 'ye(Y)*'
        else:
            y_e_string = ''

        # prefactor
        if self.prefactor != 1.0:
            prefactor_string = f"{self.prefactor:1.14e}*"
        else:
            prefactor_string = ""

        return "{}{}{}{}*rate_eval.{}".format(prefactor_string, dens_string,
                                           y_e_string, Y_string, self.fname)

    def eval(self, T, rhoY=None):
        raise NotImplementedError("base Rate class does not know how to eval()")

    def jacobian_string_py(self, y_i):
        """
        Return a string containing the term in a jacobian matrix
        in a reaction network corresponding to this rate differentiated
        with respect to y_i

        y_i is an objecs of the class ``Nucleus``.
        """
        if y_i not in self.reactants:
            return ""

        # composition dependence
        Y_string = ""
        for n, r in enumerate(sorted(set(self.reactants))):
            c = self.reactants.count(r)
            if y_i == r:
                # take the derivative
                if c == 1:
                    continue
                if 0 < n < len(set(self.reactants))-1:
                    Y_string += "*"
                if c > 2:
                    Y_string += f"{c}*Y[j{r}]**{c-1}"
                elif c == 2:
                    Y_string += f"2*Y[j{r}]"
            else:
                # this nucleus is in the rate form, but we are not
                # differentiating with respect to it
                if 0 < n < len(set(self.reactants))-1:
                    Y_string += "*"
                if c > 1:
                    Y_string += f"Y[j{r}]**{c}"
                else:
                    Y_string += f"Y[j{r}]"

        # density dependence
        if self.dens_exp == 0:
            dens_string = ""
        elif self.dens_exp == 1:
            dens_string = "rho*"
        else:
            dens_string = f"rho**{self.dens_exp}*"

        # electron fraction dependence
        if self.weak_type == 'electron_capture' and not self.tabular:
            y_e_string = 'ye(Y)*'
        else:
            y_e_string = ""

        # prefactor
        if self.prefactor != 1.0:
            prefactor_string = f"{self.prefactor:1.14e}*"
        else:
            prefactor_string = ""

        if Y_string == "" and dens_string == "" and prefactor_string == "" and y_e_string == "":
            rstring = "{}{}{}{}rate_eval.{}"
        else:
            rstring = "{}{}{}{}*rate_eval.{}"
        return rstring.format(prefactor_string, dens_string,
                              y_e_string, Y_string, self.fname)

    def eval_jacobian_term(self, T, rho, comp, y_i):
        """Evaluate drate/d(y_i), y_i is a Nucleus object.  This rate
        term has the full composition and density dependence, i.e.:

          rate = rho**n Y1**a Y2**b ... N_A <sigma v>

        The derivative is only non-zero if this term depends on
        nucleus y_i.

        """
        if y_i not in self.reactants:
            return 0.0

        ymolar = comp.get_molar()

        # composition dependence
        Y_term = 1.0
        for n, r in enumerate(sorted(set(self.reactants))):
            c = self.reactants.count(r)
            if y_i == r:
                # take the derivative
                if c == 1:
                    continue
                Y_term *= c * ymolar[r]**(c-1)
            else:
                # this nucleus is in the rate form, but we are not
                # differentiating with respect to it
                Y_term *= ymolar[r]**c

        # density dependence
        dens_term = rho**self.dens_exp

        # electron fraction dependence
        if self.weak_type == 'electron_capture' and not self.tabular:
            y_e_term = comp.eval_ye()
        else:
            y_e_term = 1.0

        # finally evaluate the rate -- for tabular rates, we need to set rhoY
        rate_eval = self.eval(T, rhoY=rho*comp.eval_ye())

        return self.prefactor * dens_term * y_e_term * Y_term * rate_eval


class ReacLibRate(Rate):
    """A single reaction rate.  Currently, this is a ReacLib rate, which
    can be composed of multiple sets, or a tabulated electron capture
    rate.

    """
    def __init__(self, rfile=None, rfile_path=None, chapter=None, original_source=None,
                 reactants=None, products=None, sets=None, labelprops=None, Q=None):
        """ rfile can be either a string specifying the path to a rate file or
        an io.StringIO object from which to read rate information. """

        self.rfile_path = rfile_path
        self.rfile = None

        if isinstance(rfile, str):
            self.rfile_path = _find_rate_file(rfile)
            self.rfile = os.path.basename(rfile)

        self.chapter = chapter    # the Reaclib chapter for this reaction
        self.original_source = original_source   # the contents of the original rate file
        self.fname = None

        if reactants:
            self.reactants = reactants
        else:
            self.reactants = []

        if products:
            self.products = products
        else:
            self.products = []

        if sets:
            self.sets = sets
        else:
            self.sets = []

        # a modified rate is one where we manually changed some of its
        # properties

        self.modified = False

        self.labelprops = labelprops

        self.approx = False
        if self.labelprops == "approx":
            self.approx = True

        self.derived = False
        if self.labelprops == "derived":
            self.derived = True

        self.label = None
        self.resonant = None
        self.weak = None
        self.weak_type = None
        self.reverse = None

        self.removed = None

        self.Q = Q

        self.tabular = False

        if isinstance(rfile, str):
            # read in the file, parse the different sets and store them as
            # SingleSet objects in sets[]
            f = open(self.rfile_path)
        elif isinstance(rfile, io.StringIO):
            # Set f to the io.StringIO object
            f = rfile
        else:
            f = None

        if f:
            self._read_from_file(f)
            f.close()
        else:
            self._set_label_properties()

        self._set_rhs_properties()
        self._set_screening()
        self._set_print_representation()

    def _set_print_representation(self):
        """ compose the string representations of this Rate. """

        super()._set_print_representation()

        if not self.fname:
            # This is used to determine which rates to detect as the same reaction
            # from multiple sources in a Library file, so it should not be unique
            # to a given source, e.g. wc12, but only unique to the reaction.
            reactants_str = '_'.join([repr(nuc) for nuc in self.reactants])
            products_str = '_'.join([repr(nuc) for nuc in self.products])
            self.fname = f'{reactants_str}__{products_str}'
            if self.weak:
                self.fname += f'__weak__{self.weak_type}'
            if self.modified:
                self.fname += "__modified"
            if self.approx:
                self.fname += "__approx"
            if self.derived:
                self.fname += "__derived"
            if self.removed:
                self.fname += "__removed"

    def modify_products(self, new_products):
        if not isinstance(new_products, (set, list, tuple)):
            new_products = [new_products]

        self.products = []
        for p in new_products:
            if isinstance(p, Nucleus):
                self.products.append(p)
            else:
                self.products.append(Nucleus(p))

        self.modified = True

        # we need to update the Q value and the print string for the rate

        self._set_q()
        self._set_screening()
        self.fname = None    # reset so it will be updated
        self._set_print_representation()

    def __hash__(self):
        return hash(self.__repr__())

    def __eq__(self, other):
        """ Determine whether two Rate objects are equal.
        They are equal if they contain identical reactants and products and
        if they contain the same SingleSet sets and if their chapters are equal."""

        if not isinstance(other, ReacLibRate):
            return False

        x = True

        x = x and (self.chapter == other.chapter)
        x = x and (self.reactants == other.reactants)
        x = x and (self.products == other.products)
        x = x and (len(self.sets) == len(other.sets))

        for si in self.sets:
            scomp = False
            for sj in other.sets:
                if si == sj:
                    scomp = True
                    break
            x = x and scomp

        return x

    def __add__(self, other):
        """Combine the sets of two Rate objects if they describe the same
           reaction. Must be Reaclib rates."""
        assert self.reactants == other.reactants
        assert self.products == other.products
        assert self.chapter == other.chapter
        assert isinstance(self.chapter, int)
        assert self.label == other.label
        assert self.weak == other.weak
        assert self.weak_type == other.weak_type
        assert self.reverse == other.reverse

        if self.resonant != other.resonant:
            self._labelprops_combine_resonance()
        new_rate = ReacLibRate(chapter=self.chapter,
                               original_source='\n'.join([self.original_source,
                                                          other.original_source]),
                               reactants=self.reactants,
                               products=self.products,
                               sets=self.sets + other.sets,
                               labelprops=self.labelprops,
                               Q=self.Q)
        return new_rate

    def _set_label_properties(self, labelprops=None):
        """ Calls _update_resonance_combined and then
            _update_label_properties. """
        if labelprops:
            self.labelprops = labelprops

        # Update labelprops based on the Sets in this Rate
        # to set the resonance_combined flag properly
        self._update_resonance_combined()
        self._update_label_properties()

    def _update_resonance_combined(self):
        """ Checks the Sets in this Rate and updates the
            resonance_combined flag as well as
            self.labelprops[4] """
        sres = [s.resonant for s in self.sets]
        if True in sres and False in sres:
            self._labelprops_combine_resonance()

    def _labelprops_combine_resonance(self):
        """ Update self.labelprops[4] = 'c'"""
        llp = list(self.labelprops)
        llp[4] = 'c'
        self.labelprops = ''.join(llp)

    def _update_label_properties(self):
        """ Set label and flags indicating Rate is resonant,
            weak, or reverse. """
        assert isinstance(self.labelprops, str)
        if self.labelprops == "approx":
            self.label = "approx"
            self.resonant = False
            self.weak = False
            self.weak_type = None
            self.reverse = False
        elif self.labelprops == "derived":
            self.label = "derived"
            self.resonant = False  # Derived may be resonant in some cases
            self.weak = False
            self.weak_type = None
            self.reverse = False
        else:
            assert len(self.labelprops) == 6
            self.label = self.labelprops[0:4]
            self.resonant = self.labelprops[4] == 'r'
            self.weak = self.labelprops[4] == 'w'
            if self.weak:
                if self.label.strip() == 'ec' or self.label.strip() == 'bec':
                    self.weak_type = 'electron_capture'
                else:
                    self.weak_type = self.label.strip().replace('+', '_pos_').replace('-', '_neg_')
            else:
                self.weak_type = None
            self.reverse = self.labelprops[5] == 'v'

    def _read_from_file(self, f):
        """ given a file object, read rate data from the file. """
        lines = f.readlines()
        f.close()

        self.original_source = "".join(lines)

        # first line is the chapter
        self.chapter = lines[0].strip()
        self.chapter = int(self.chapter)

        # remove any blank lines
        set_lines = [l for l in lines[1:] if not l.strip() == ""]

        # the rest is the sets
        first = 1
        while len(set_lines) > 0:
            # check for a new chapter id in case of Reaclib v2 format
            check_chapter = set_lines[0].strip()
            try:
                # see if there is a chapter number preceding the set
                check_chapter = int(check_chapter)
                # check that the chapter number is the same as the first
                # set in this rate file
                assert check_chapter == self.chapter, f'read chapter {check_chapter}, expected chapter {self.chapter} for this rate set.'
                # get rid of chapter number so we can read a rate set
                set_lines.pop(0)
            except (TypeError, ValueError):
                # there was no chapter number, proceed reading a set
                pass

            # sets are 3 lines long
            s1 = set_lines.pop(0)
            s2 = set_lines.pop(0)
            s3 = set_lines.pop(0)

            # first line of a set has up to 6 nuclei, then the label,
            # and finally the Q value

            # get rid of first 5 spaces
            s1 = s1[5:]

            # next follows 6 fields of 5 characters containing nuclei
            # the 6 fields are padded with spaces
            f = []
            for i in range(6):
                ni = s1[:5]
                s1 = s1[5:]
                ni = ni.strip()
                if ni:
                    f.append(ni)

            # next come 8 spaces, so get rid of them
            s1 = s1[8:]

            # next is a 4-character set label and 2 character flags
            labelprops = s1[:6]
            s1 = s1[6:]

            # next come 3 spaces
            s1 = s1[3:]

            # next comes a 12 character Q value followed by 10 spaces
            Q = float(s1.strip())

            if first:
                self.Q = Q

                # what's left are the nuclei -- their interpretation
                # depends on the chapter
                if self.chapter == 1:
                    # e1 -> e2
                    self.reactants.append(Nucleus.from_cache(f[0]))
                    self.products.append(Nucleus.from_cache(f[1]))

                elif self.chapter == 2:
                    # e1 -> e2 + e3
                    self.reactants.append(Nucleus.from_cache(f[0]))
                    self.products += [Nucleus.from_cache(f[1]), Nucleus.from_cache(f[2])]

                elif self.chapter == 3:
                    # e1 -> e2 + e3 + e4
                    self.reactants.append(Nucleus.from_cache(f[0]))
                    self.products += [Nucleus.from_cache(f[1]), Nucleus.from_cache(f[2]),
                                      Nucleus.from_cache(f[3])]

                elif self.chapter == 4:
                    # e1 + e2 -> e3
                    self.reactants += [Nucleus.from_cache(f[0]), Nucleus.from_cache(f[1])]
                    self.products.append(Nucleus.from_cache(f[2]))

                elif self.chapter == 5:
                    # e1 + e2 -> e3 + e4
                    self.reactants += [Nucleus.from_cache(f[0]), Nucleus.from_cache(f[1])]
                    self.products += [Nucleus.from_cache(f[2]), Nucleus.from_cache(f[3])]

                elif self.chapter == 6:
                    # e1 + e2 -> e3 + e4 + e5
                    self.reactants += [Nucleus.from_cache(f[0]), Nucleus.from_cache(f[1])]
                    self.products += [Nucleus.from_cache(f[2]), Nucleus.from_cache(f[3]),
                                      Nucleus.from_cache(f[4])]

                elif self.chapter == 7:
                    # e1 + e2 -> e3 + e4 + e5 + e6
                    self.reactants += [Nucleus.from_cache(f[0]), Nucleus.from_cache(f[1])]
                    self.products += [Nucleus.from_cache(f[2]), Nucleus.from_cache(f[3]),
                                      Nucleus.from_cache(f[4]), Nucleus.from_cache(f[5])]

                elif self.chapter == 8:
                    # e1 + e2 + e3 -> e4
                    self.reactants += [Nucleus.from_cache(f[0]), Nucleus.from_cache(f[1]),
                                       Nucleus.from_cache(f[2])]
                    self.products.append(Nucleus.from_cache(f[3]))
                    # support historical format, where chapter 8 also handles what are
                    # now chapter 9 rates
                    if len(f) == 5:
                        self.products.append(Nucleus.from_cache(f[4]))

                elif self.chapter == 9:
                    # e1 + e2 + e3 -> e4 + e5
                    self.reactants += [Nucleus.from_cache(f[0]), Nucleus.from_cache(f[1]),
                                       Nucleus.from_cache(f[2])]
                    self.products += [Nucleus.from_cache(f[3]), Nucleus.from_cache(f[4])]

                elif self.chapter == 10:
                    # e1 + e2 + e3 + e4 -> e5 + e6
                    self.reactants += [Nucleus.from_cache(f[0]), Nucleus.from_cache(f[1]),
                                       Nucleus.from_cache(f[2]), Nucleus.from_cache(f[3])]
                    self.products += [Nucleus.from_cache(f[4]), Nucleus.from_cache(f[5])]

                elif self.chapter == 11:
                    # e1 -> e2 + e3 + e4 + e5
                    self.reactants.append(Nucleus.from_cache(f[0]))
                    self.products += [Nucleus.from_cache(f[1]), Nucleus.from_cache(f[2]),
                                      Nucleus.from_cache(f[3]), Nucleus.from_cache(f[4])]
                else:
                    print(f'Chapter could not be identified in {self.original_source}')
                    assert isinstance(self.chapter, int) and self.chapter <= 11

                first = 0

            # the second line contains the first 4 coefficients
            # the third lines contains the final 3
            # we can't just use split() here, since the fields run into one another
            n = 13  # length of the field
            a = [s2[i:i+n] for i in range(0, len(s2), n)]
            a += [s3[i:i+n] for i in range(0, len(s3), n)]

            a = [float(e) for e in a if not e.strip() == ""]
            self.sets.append(SingleSet(a, labelprops=labelprops))
            self._set_label_properties(labelprops)

    def get_rate_id(self):
        """ Get an identifying string for this rate.
        Don't include resonance state since we combine resonant and
        non-resonant versions of reactions. """

        srev = ''
        if self.reverse:
            srev = 'reverse'

        sweak = ''
        if self.weak:
            sweak = 'weak'

        ssrc = 'reaclib'

        return f'{self.rid} <{self.label.strip()}_{ssrc}_{sweak}_{srev}>'

    def function_string_py(self):
        """
        Return a string containing python function that computes the
        rate
        """

        fstring = ""
        fstring += "@numba.njit()\n"
        fstring += f"def {self.fname}(rate_eval, tf):\n"
        fstring += f"    # {self.rid}\n"
        fstring += "    rate = 0.0\n\n"

        for s in self.sets:
            fstring += f"    # {s.labelprops[0:5]}\n"
            set_string = s.set_string_py(prefix="rate", plus_equal=True)
            for t in set_string.split("\n"):
                fstring += "    " + t + "\n"

        fstring += "\n"
        fstring += f"    rate_eval.{self.fname} = rate\n\n"
        return fstring

    def function_string_cxx(self, dtype="double", specifiers="inline", leave_open=False):
        """
        Return a string containing C++ function that computes the
        rate
        """

        fstring = ""
        fstring += "template <int do_T_derivatives>\n"
        fstring += f"{specifiers}\n"
        fstring += f"void rate_{self.fname}(const tf_t& tfactors, {dtype}& rate, {dtype}& drate_dT) {{\n\n"
        fstring += f"    // {self.rid}\n\n"
        fstring += "    rate = 0.0;\n"
        fstring += "    drate_dT = 0.0;\n\n"
        fstring += f"    {dtype} ln_set_rate{{0.0}};\n"
        fstring += f"    {dtype} dln_set_rate_dT9{{0.0}};\n"
        fstring += f"    {dtype} set_rate{{0.0}};\n\n"

        for s in self.sets:
            fstring += f"    // {s.labelprops[0:5]}\n"
            set_string = s.set_string_cxx(prefix="ln_set_rate", plus_equal=False, with_exp=False)
            for t in set_string.split("\n"):
                fstring += "    " + t + "\n"
            fstring += "\n"

            fstring += "    if constexpr (do_T_derivatives) {\n"
            dln_set_string_dT9 = s.dln_set_string_dT9_cxx(prefix="dln_set_rate_dT9", plus_equal=False)
            for t in dln_set_string_dT9.split("\n"):
                fstring += "        " + t + "\n"
            fstring += "    }\n"
            fstring += "\n"

            fstring += "    // avoid underflows by zeroing rates in [0.0, 1.e-100]\n"
            fstring += "    ln_set_rate = std::max(ln_set_rate, -230.0);\n"
            fstring += "    set_rate = std::exp(ln_set_rate);\n"

            fstring += "    rate += set_rate;\n"

            fstring += "    if constexpr (do_T_derivatives) {\n"
            fstring += "        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;\n"
            fstring += "    }\n\n"

        if not leave_open:
            fstring += "}\n\n"

        return fstring

    def eval(self, T, rhoY=None):
        """ evauate the reaction rate for temperature T """

        tf = Tfactors(T)
        r = 0.0
        for s in self.sets:
            f = s.f()
            r += f(tf)

        return r

    def eval_deriv(self, T, rhoY=None):
        """ evauate the derivative of reaction rate with respect to T """

        tf = Tfactors(T)
        drdT = 0.0
        for s in self.sets:
            dfdT = s.dfdT()
            drdT += dfdT(tf)

        return drdT

    def get_rate_exponent(self, T0):
        """
        for a rate written as a power law, r = r_0 (T/T0)**nu, return
        nu corresponding to T0
        """

        # nu = dln r /dln T, so we need dr/dT
        r1 = self.eval(T0)
        dT = 1.e-8*T0
        r2 = self.eval(T0 + dT)

        drdT = (r2 - r1)/dT
        return (T0/r1)*drdT

    def plot(self, Tmin=1.e8, Tmax=1.6e9, rhoYmin=3.9e8, rhoYmax=2.e9,
             figsize=(10, 10)):
        """plot the rate's temperature sensitivity vs temperature

        :param float Tmin:    minimum temperature for plot
        :param float Tmax:    maximum temperature for plot
        :param float rhoYmin: minimum electron density to plot (e-capture rates only)
        :param float rhoYmax: maximum electron density to plot (e-capture rates only)
        :param tuple figsize: figure size specification for matplotlib

        :return: a matplotlib figure object
        :rtype: matplotlib.figure.Figure

        """

        fig, ax = plt.subplots(figsize=figsize)

        temps = np.logspace(np.log10(Tmin), np.log10(Tmax), 100)
        r = np.zeros_like(temps)

        for n, T in enumerate(temps):
            r[n] = self.eval(T)

        ax.loglog(temps, r)
        ax.set_xlabel(r"$T$")

        if self.dens_exp == 0:
            ax.set_ylabel(r"\tau")
        elif self.dens_exp == 1:
            ax.set_ylabel(r"$N_A <\sigma v>$")
        elif self.dens_exp == 2:
            ax.set_ylabel(r"$N_A^2 <n_a n_b n_c v>$")

        ax.set_title(fr"{self.pretty_string}")

        return fig


class TabularRate(Rate):
    """A tabular rate."""
    def __init__(self, rfile=None, rfile_path=None):
        """ rfile can be either a string specifying the path to a rate file or
        an io.StringIO object from which to read rate information. """

        self.rfile_path = rfile_path
        self.rfile = None

        if isinstance(rfile, str):
            self.rfile_path = _find_rate_file(rfile)
            self.rfile = os.path.basename(rfile)

        self.fname = None

        self.reactants = []
        self.products = []

        self.label = "tabular"
        self.tabular = True

        self.Q = None

        # we should initialize this somehow
        self.weak_type = ""

        if isinstance(rfile, str):
            # read in the file, parse the different sets and store them as
            # SingleSet objects in sets[]
            f = open(self.rfile_path)
        elif isinstance(rfile, io.StringIO):
            # Set f to the io.StringIO object
            f = rfile
        else:
            f = None

        if f:
            self._read_from_file(f)
            f.close()

        self._set_rhs_properties()
        self._set_screening()
        self._set_print_representation()

        self.get_tabular_rate()

    def __hash__(self):
        return hash(self.__repr__())

    def __eq__(self, other):
        """ Determine whether two Rate objects are equal.
        They are equal if they contain identical reactants and products."""

        if not isinstance(other, TabularRate):
            return False

        return self.reactants == other.reactants and self.products == other.products

    def __add__(self, other):
        raise NotImplementedError("addition not defined for tabular rates")

    def _read_from_file(self, f):
        """ given a file object, read rate data from the file. """
        lines = f.readlines()
        f.close()

        self.original_source = "".join(lines)

        # first line is the chapter
        self.chapter = lines[0].strip()
        assert self.chapter == "t"

        # remove any blank lines
        set_lines = [l for l in lines[1:] if not l.strip() == ""]

        # e1 -> e2, Tabulated
        s1 = set_lines.pop(0)
        s2 = set_lines.pop(0)
        s3 = set_lines.pop(0)
        s4 = set_lines.pop(0)
        s5 = set_lines.pop(0)
        f = s1.split()
        try:
            self.reactants.append(Nucleus.from_cache(f[0]))
            self.products.append(Nucleus.from_cache(f[1]))
        except Exception as ex:
            raise Exception(f'Nucleus objects could not be identified in {self.original_source}') from ex

        self.table_file = s2.strip()
        self.table_header_lines = int(s3.strip())
        self.table_rhoy_lines = int(s4.strip())
        self.table_temp_lines = int(s5.strip())
        self.table_num_vars = 6  # Hard-coded number of variables in tables for now.
        self.table_index_name = f'j_{self.reactants[0]}_{self.products[0]}'
        self.labelprops = 'tabular'

        # set weak type
        if "electroncapture" in self.table_file:
            self.weak_type = "electron_capture"

        elif "betadecay" in self.table_file:
            self.weak_type = "beta_decay"

    def _set_rhs_properties(self):
        """ compute statistical prefactor and density exponent from the reactants. """
        self.prefactor = 1.0  # this is 1/2 for rates like a + a (double counting)
        self.inv_prefactor = 1
        for r in set(self.reactants):
            self.inv_prefactor = self.inv_prefactor * np.math.factorial(self.reactants.count(r))
        self.prefactor = self.prefactor/float(self.inv_prefactor)
        self.dens_exp = len(self.reactants)-1

    def _set_screening(self):
        """ tabular rates are not currently screened (they are e-capture or beta-decay)"""
        self.ion_screen = []
        self.symmetric_screen = []

        if not self.fname:
            # This is used to determine which rates to detect as the same reaction
            # from multiple sources in a Library file, so it should not be unique
            # to a given source, e.g. wc12, but only unique to the reaction.
            reactants_str = '_'.join([repr(nuc) for nuc in self.reactants])
            products_str = '_'.join([repr(nuc) for nuc in self.products])
            self.fname = f'{reactants_str}__{products_str}'

    def get_rate_id(self):
        """ Get an identifying string for this rate.
        Don't include resonance state since we combine resonant and
        non-resonant versions of reactions. """

        ssrc = 'tabular'

        return f'{self.rid} <{self.label.strip()}_{ssrc}>'

    def function_string_py(self):
        """
        Return a string containing python function that computes the
        rate
        """

        fstring = ""
        fstring += "@numba.njit()\n"
        fstring += f"def {self.fname}(rate_eval, T, rhoY):\n"
        fstring += f"    # {self.rid}\n"

        # find the nearest value of T and rhoY in the data table
        fstring += f"    T_nearest = ({self.fname}_data[:, 1])[np.abs(({self.fname}_data[:, 1]) - T).argmin()]\n"
        fstring += f"    rhoY_nearest = ({self.fname}_data[:, 0])[np.abs(({self.fname}_data[:, 0]) - rhoY).argmin()]\n"
        fstring += f"    inde = np.where(({self.fname}_data[:, 1] == T_nearest) & ({self.fname}_data[:, 0] == rhoY_nearest))[0][0]\n"
        fstring += f"    rate_eval.{self.fname} = {self.fname}_data[inde][5]\n\n"

        return fstring

    def get_tabular_rate(self):
        """read the rate data from .dat file """

        # find .dat file and read it
        self.table_path = _find_rate_file(self.table_file)
        t_data2d = []
        with open(self.table_path) as tabular_file:
            for i, line in enumerate(tabular_file):
                # skip header lines
                if i < self.table_header_lines:
                    continue
                line = line.strip()
                # skip empty lines
                if not line:
                    continue
                # split the column values on whitespace
                t_data2d.append(line.split())

        # convert the nested list of string values into a numpy float array
        self.tabular_data_table = np.array(t_data2d, dtype=float)

    def eval(self, T, rhoY=None):
        """ evauate the reaction rate for temperature T """

        data = self.tabular_data_table
        # find the nearest value of T and rhoY in the data table
        T_nearest = (data[:, 1])[np.abs((data[:, 1]) - T).argmin()]
        rhoY_nearest = (data[:, 0])[np.abs((data[:, 0]) - rhoY).argmin()]
        inde = np.where((data[:, 1] == T_nearest) & (data[:, 0] == rhoY_nearest))[0][0]
        r = data[inde][5]
        return r

    def get_nu_loss(self, T, rhoY):
        """ get the neutrino loss rate for the reaction if tabulated"""

        nu_loss = None
        data = self.tabular_data_table
        # find the nearest value of T and rhoY in the data table
        T_nearest = (data[:, 1])[np.abs((data[:, 1]) - T).argmin()]
        rhoY_nearest = (data[:, 0])[np.abs((data[:, 0]) - rhoY).argmin()]
        inde = np.where((data[:, 1] == T_nearest) & (data[:, 0] == rhoY_nearest))[0][0]
        nu_loss = data[inde][6]

        return nu_loss

    def plot(self, Tmin=1.e8, Tmax=1.6e9, rhoYmin=3.9e8, rhoYmax=2.e9, color_field='rate',
             figsize=(10, 10)):
        """plot the rate's temperature sensitivity vs temperature

        :param float Tmin:    minimum temperature for plot
        :param float Tmax:    maximum temperature for plot
        :param float rhoYmin: minimum electron density to plot (e-capture rates only)
        :param float rhoYmax: maximum electron density to plot (e-capture rates only)
        :param tuple figsize: figure size specification for matplotlib

        :return: a matplotlib figure object
        :rtype: matplotlib.figure.Figure

        """

        fig, ax = plt.subplots(figsize=figsize)

        data = self.tabular_data_table

        inde1 = data[:, 1] <= Tmax
        inde2 = data[:, 1] >= Tmin
        inde3 = data[:, 0] <= rhoYmax
        inde4 = data[:, 0] >= rhoYmin
        data_heatmap = data[inde1 & inde2 & inde3 & inde4].copy()

        rows, row_pos = np.unique(data_heatmap[:, 0], return_inverse=True)
        cols, col_pos = np.unique(data_heatmap[:, 1], return_inverse=True)
        pivot_table = np.zeros((len(rows), len(cols)), dtype=data_heatmap.dtype)

        if color_field == 'rate':
            icol = 5
            title = f"{self.weak_type} rate in log10(1/s)"
            cmap = 'magma'

        elif color_field == 'nu_loss':
            icol = 6
            title = "neutrino energy loss rate in log10(erg/s)"
            cmap = 'viridis'

        else:
            raise ValueError("color_field must be either 'rate' or 'nu_loss'.")

        try:
            pivot_table[row_pos, col_pos] = np.log10(data_heatmap[:, icol])
        except ValueError:
            print("Divide by zero encountered in log10\nChange the scale of T or rhoY")

        im = ax.imshow(pivot_table, cmap=cmap)
        fig.colorbar(im, ax=ax)

        ax.set_xlabel(r"$\log(T)$ [K]")
        ax.set_ylabel(r"$\log(\rho Y_e)$ [g/cm$^3$]")
        ax.set_title(fr"{self.pretty_string}" + "\n" + title)
        ax.set_yticks(range(len(rows)))
        ylabels = [f"{np.log10(q):4.2f}" for q in rows]
        ax.set_yticklabels(ylabels)
        ax.set_xticks(range(len(cols)))
        xlabels = [f"{np.log10(q):4.2f}" for q in cols]
        ax.set_xticklabels(xlabels, rotation=90, ha="right", rotation_mode="anchor")
        ax.invert_yaxis()

        return fig


class DerivedRate(ReacLibRate):
    """
    This class is a derived class from `Rate` with the purpose of computing the inverse rate
    by the application of detailed balance to the forward reactions.
    """

    def __init__(self, rate,  compute_Q=False, use_pf=False):

        self.use_pf = use_pf
        self.rate = rate
        self.compute_Q = compute_Q

        assert isinstance(rate, Rate)

        if (isinstance(rate, TabularRate) or self.rate.weak or
            self.rate.reverse):
            raise ValueError('The rate is reverse or weak or tabular')

        for nuc in self.rate.reactants:

            if not nuc.spin_states:
                raise ValueError('One of the reactants spin ground state, is not defined')

        for nuc in self.rate.products:

            if not nuc.spin_states:
                raise ValueError('One of the products spin ground state, is not defined')

        derived_sets = []
        for ssets in self.rate.sets:
            a = ssets.a
            prefactor = 0.0
            Q = 0.0
            prefactor += -np.log(N_a) * (len(self.rate.reactants) - len(self.rate.products))

            for nucr in self.rate.reactants:
                prefactor += 1.5*np.log(nucr.A) + np.log(nucr.spin_states)
                Q += nucr.A_nuc
            for nucp in self.rate.products:
                prefactor += -1.5*np.log(nucp.A) - np.log(nucp.spin_states)
                Q -= nucp.A_nuc

            if self.compute_Q:
                Q = Q * amu_mev
            else:
                Q = self.rate.Q

            prefactor += np.log(self.counter_factors()[1]) - np.log(self.counter_factors()[0])

            if len(self.rate.reactants) == len(self.rate.products):
                prefactor += 0.0
            else:
                F = (amu * k_B * 1.0e5 / (2.0*np.pi*hbar**2))**(1.5*(len(self.rate.reactants) - len(self.rate.products)))
                prefactor += np.log(F)

            a_rev = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            a_rev[0] = prefactor + a[0]
            a_rev[1] = a[1] - Q / (1.0e9 * k_B_mev_k)
            a_rev[2] = a[2]
            a_rev[3] = a[3]
            a_rev[4] = a[4]
            a_rev[5] = a[5]
            a_rev[6] = a[6] + 1.5*(len(self.rate.reactants) - len(self.rate.products))
            sset_d = SingleSet(a=a_rev, labelprops=rate.labelprops)
            derived_sets.append(sset_d)

        super().__init__(rfile=self.rate.rfile, rfile_path=self.rate.rfile_path, chapter=self.rate.chapter, original_source=self.rate.original_source,
                reactants=self.rate.products, products=self.rate.reactants, sets=derived_sets, labelprops="derived", Q=-Q)

    def eval(self, T, rhoY=None):

        r = super().eval(T=T, rhoY=rhoY)
        z_r = 1.0
        z_p = 1.0
        if self.use_pf:
            for nuc in set(self.rate.reactants + self.rate.products):
                if not nuc.partition_function and str(nuc) != 'h1' and str(nuc) != 'n' and str(nuc) != 'he4' and str(nuc) != 'p':
                    print(f'WARNING: {nuc} partition function is not supported by tables: set pf = 1.0 by default')

            for nucr in self.rate.reactants:
                if not nucr.partition_function:
                    continue
                    #nucr.partition_function = lambda T: 1.0
                z_r *= nucr.partition_function(T)

            for nucp in self.rate.products:
                if not nucp.partition_function:
                    continue
                    #nucp.partition_function = lambda T: 1.0
                z_p *= nucp.partition_function(T)

            return r*z_r/z_p
        return r

    def function_string_py(self):
        """
        Return a string containing python function that computes the
        rate
        """

        for nuc in set(self.rate.reactants + self.rate.products):
            if not nuc.partition_function and str(nuc) not in ['h1', 'n', 'he4', 'p']:
                print(f'WARNING: {nuc} partition function is not supported by tables: set pf = 1.0 by default')

        fstring = super().function_string_py()

        if self.use_pf:

            fstring += "\n"
            for nuc in set(self.rate.reactants + self.rate.products):
                if nuc.partition_function:
                    fstring += f"    # interpolating {nuc} partition function\n"
                    fstring += f"    {nuc}_pf_exponent = np.interp(tf.T9, xp={nuc}_temp_array, fp=np.log10({nuc}_pf_array))\n"
                    fstring += f"    {nuc}_pf = 10.0**{nuc}_pf_exponent\n"
                else:
                    fstring += f"    # setting {nuc} partition function to 1.0 by default, independent of T\n"
                    fstring += f"    {nuc}_pf = 1.0\n"
                fstring += "\n"

            fstring += "    "
            fstring += "z_r = "
            fstring += "*".join([f"{nucr}_pf" for nucr in self.rate.reactants])

            fstring += "\n"
            fstring += "    "
            fstring += "z_p = "
            fstring += "*".join([f"{nucp}_pf" for nucp in self.rate.products])

            fstring += "\n"
            fstring += f"    rate_eval.{self.fname} *= z_r/z_p\n"

        return fstring

    def function_string_cxx(self, dtype="double", specifiers="inline", leave_open=False):
        """
        Return a string containing C++ function that computes the
        rate
        """

        for nuc in set(self.rate.reactants + self.rate.products):
            if not nuc.partition_function and str(nuc) not in ['h1', 'n', 'he4', 'p']:
                print(f'WARNING: {nuc} partition function is not supported by tables: set pf = 1.0 by default')

        fstring = super().function_string_cxx(dtype=dtype, specifiers=specifiers, leave_open=True)

        # right now we have rate and drate_dT without the partition function
        # now the partition function corrections

        if self.use_pf:

            fstring += "\n"
            for nuc in set(self.rate.reactants + self.rate.products):
                fstring += f"    Real {nuc}_pf, d{nuc}_pf_dT;\n"

                if nuc.partition_function:
                    fstring += f"    // interpolating {nuc} partition function\n"
                    fstring += f"    get_partition_function({nuc.cindex()}, tfactors, {nuc}_pf, d{nuc}_pf_dT);\n"
                else:
                    fstring += f"    // setting {nuc} partition function to 1.0 by default, independent of T\n"
                    fstring += f"    {nuc}_pf = 1.0_rt;\n"
                    fstring += f"    d{nuc}_pf_dT = 0.0_rt;\n"
                fstring += "\n"

            fstring += "    Real z_r = "
            fstring += " * ".join([f"{nucr}_pf" for nucr in self.rate.reactants])
            fstring += ";\n"

            fstring += "    Real z_p = "
            fstring += " * ".join([f"{nucp}_pf" for nucp in self.rate.products])
            fstring += ";\n\n"

            # now the derivatives, via chain rule
            chain_terms = []
            for n in self.rate.reactants:
                chain_terms.append(" * ".join([f"{nucr}_pf" for nucr in self.rate.reactants if nucr != n] + [f"d{n}_pf_dT"]))

            fstring += "    Real dz_r_dT = "
            fstring += " + ".join(chain_terms)
            fstring += ";\n"

            chain_terms = []
            for n in self.rate.products:
                chain_terms.append(" * ".join([f"{nucp}_pf" for nucp in self.rate.products if nucp != n] + [f"d{n}_pf_dT"]))

            fstring += "    Real dz_p_dT = "
            fstring += " + ".join(chain_terms)
            fstring += ";\n\n"

            fstring += "    Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);\n\n"

            # final terms

            fstring += "    drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);\n"
            fstring += "    rate *= z_r/z_p;\n\n"

        if not leave_open:
            fstring += "}\n\n"

        return fstring

    def counter_factors(self):
        """This function returns the nucr! = nucr_1! * ... * nucr_r!
        for each repeated nucr reactant and nucp! = nucp_1! * ... *
        nucp_p! for each reactant nucp product in a ordered pair
        (nucr!, nucp!). The factors nucr! and nucp! avoid overcounting
        when more than one nuclei is involve in the reaction,
        otherwise it will return 1.0.

        """

        react_counts = Counter(self.rate.reactants)
        prod_counts = Counter(self.rate.products)

        reactant_factor = 1.0
        for nuc in set(self.rate.reactants):
            reactant_factor *= np.math.factorial(react_counts[nuc])

        product_factor = 1.0
        for nuc in set(self.rate.products):
            product_factor *= np.math.factorial(prod_counts[nuc])

        return (reactant_factor, product_factor)


class RatePair:
    """the forward and reverse rates for a single reaction sequence.
    Forward rates are those with Q >= 0.

    :var forward: the forward reaction Rate object
    :var reverse: the reverse reaction Rate object

    """

    def __init__(self, forward=None, reverse=None):
        self.forward = forward
        self.reverse = reverse

    def __repr__(self):
        return f"forward: {self.forward} ; reverse: {self.reverse}"

    def __lt__(self, other):
        if self.forward is not None and other.forward is not None:
            return self.forward < other.forward
        if self.forward is None:
            return False
        return True

    def __eq__(self, other):
        return self.forward == other.forward and self.reverse == other.reverse


class ApproximateRate(ReacLibRate):

    def __init__(self, primary_rate, secondary_rates,
                 primary_reverse, secondary_reverse, is_reverse=False, approx_type="ap_pg"):
        """the primary rate has the same reactants and products and the final
        approximate rate would have.  The secondary rates are ordered such that
        together they would give the same sequence"""

        self.primary_rate = primary_rate
        self.secondary_rates = secondary_rates

        self.primary_reverse = primary_reverse
        self.secondary_reverse = secondary_reverse

        self.is_reverse = is_reverse

        self.approx_type = approx_type

        if self.approx_type == "ap_pg":

            # an ap_pg approximate rate combines A(a,g)B and A(a,p)X(p,g)B into a
            # single effective rate by assuming proton equilibrium.

            assert len(secondary_rates) == 2

            # make sure that the primary forward rate makes sense
            # this should be A(a,g)B

            assert Nucleus("he4") in self.primary_rate.reactants and len(self.primary_rate.products) == 1

            # we are going to define the product A and reactant B from this reaction

            self.primary_reactant = max(self.primary_rate.reactants)
            self.primary_product = max(self.primary_rate.products)

            # the first secondary rate should be A(a,p)X, where X is the
            # intermediate nucleus

            assert (self.primary_reactant in self.secondary_rates[0].reactants and
                    Nucleus("he4") in self.secondary_rates[0].reactants and
                    Nucleus("p") in self.secondary_rates[0].products)

            # the intermediate nucleus is not in our network, so make it
            # dummy

            self.intermediate_nucleus = max(self.secondary_rates[0].products)
            #self.intermediate_nucleus.dummy = True

            # now the second secondary rate show be X(p,g)B

            assert (self.intermediate_nucleus in self.secondary_rates[1].reactants and
                    Nucleus("p") in self.secondary_rates[1].reactants and
                    self.primary_product in secondary_rates[1].products)

            # now ensure that the reverse rate makes sense

            # the primary reverse rate is B(g,a)A

            assert (self.primary_product in self.primary_reverse.reactants and
                    self.primary_reactant in self.primary_reverse.products)

            # now the first secondary reverse rate should be B(g,p)X

            assert (self.primary_product in self.secondary_reverse[0].reactants and
                    self.intermediate_nucleus in secondary_reverse[0].products and
                    Nucleus("p") in secondary_reverse[0].products)

            # and the second secondary reverse rate should be X(p,a)A

            assert (self.intermediate_nucleus in self.secondary_reverse[1].reactants and
                    Nucleus("p") in self.secondary_reverse[1].reactants and
                    self.primary_reactant in self.secondary_reverse[1].products and
                    Nucleus("he4") in self.secondary_reverse[1].products)

            # now initialize the super class with these reactants and products

            if not self.is_reverse:
                super().__init__(reactants=[self.primary_reactant, Nucleus("he4")],
                                 products=[self.primary_product],
                                 labelprops="approx", chapter=-1)
            else:
                super().__init__(reactants=[self.primary_product],
                                 products=[self.primary_reactant, Nucleus("he4")],
                                 labelprops="approx", chapter=-1)

        else:
            raise NotImplementedError(f"approximation type {self.approx_type} not supported")

        # update the Q value
        self._set_q()

    def get_child_rates(self):
        """return a list of all of the rates that are used in this approximation"""
        tlist = [self.primary_rate]
        tlist += self.secondary_rates
        tlist += [self.primary_reverse]
        tlist += self.secondary_reverse
        return tlist

    def _set_screening(self):
        # the individual rates are screened -- we don't screen the combination of them
        pass

    def eval(self, T, rhoY=None):
        """evaluate the approximate rate"""

        if self.approx_type == "ap_pg":
            if not self.is_reverse:  # pylint: disable=no-else-return
                # the approximate forward rate is r_ag + r_ap r_pg / (r_pg + r_pa)
                r_ag = self.primary_rate.eval(T)
                r_ap = self.secondary_rates[0].eval(T)
                r_pg = self.secondary_rates[1].eval(T)

                r_pa = self.secondary_reverse[1].eval(T)

                return r_ag + r_ap * r_pg / (r_pg + r_pa)

            else:
                # the approximate reverse rate is r_ga + r_pa r_gp / (r_pg + r_pa)

                r_ga = self.primary_reverse.eval(T)
                r_gp = self.secondary_reverse[0].eval(T)
                r_pa = self.secondary_reverse[1].eval(T)

                r_pg = self.secondary_rates[1].eval(T)

                return r_ga + r_pa * r_gp / (r_pg + r_pa)
        raise NotImplementedError(f"approximation type {self.approx_type} not supported")

    def function_string_py(self):
        """
        Return a string containing python function that computes the
        approximate rate
        """

        if self.approx_type != "ap_pg":
            raise NotImplementedError("don't know how to work with this approximation")

        string = ""
        string += "@numba.njit()\n"
        string += f"def {self.fname}(rate_eval, tf):\n"

        if not self.is_reverse:

            # first we need to get all of the rates that make this up
            string += f"    r_ag = rate_eval.{self.primary_rate.fname}\n"
            string += f"    r_ap = rate_eval.{self.secondary_rates[0].fname}\n"
            string += f"    r_pg = rate_eval.{self.secondary_rates[1].fname}\n"
            string += f"    r_pa = rate_eval.{self.secondary_reverse[1].fname}\n"

            # now the approximation
            string += "    rate = r_ag + r_ap * r_pg / (r_pg + r_pa)\n"

        else:

            # first we need to get all of the rates that make this up
            string += f"    r_ga = rate_eval.{self.primary_reverse.fname}\n"
            string += f"    r_pa = rate_eval.{self.secondary_reverse[1].fname}\n"
            string += f"    r_gp = rate_eval.{self.secondary_reverse[0].fname}\n"
            string += f"    r_pg = rate_eval.{self.secondary_rates[1].fname}\n"

            # now the approximation
            string += "    rate = r_ga + r_pa * r_gp / (r_pg + r_pa)\n"

        string += f"    rate_eval.{self.fname} = rate\n\n"
        return string

    def function_string_cxx(self, dtype="double", specifiers="inline", leave_open=False):
        """
        Return a string containing C++ function that computes the
        approximate rate
        """

        if self.approx_type != "ap_pg":
            raise NotImplementedError("don't know how to work with this approximation")

        fstring = ""
        fstring = "template <typename T>\n"
        fstring += f"{specifiers}\n"
        fstring += f"void rate_{self.fname}(const T& rate_eval, {dtype}& rate, {dtype}& drate_dT) {{\n\n"

        if not self.is_reverse:

            # first we need to get all of the rates that make this up
            fstring += f"    {dtype} r_ag = rate_eval.screened_rates(k_{self.primary_rate.fname});\n"
            fstring += f"    {dtype} r_ap = rate_eval.screened_rates(k_{self.secondary_rates[0].fname});\n"
            fstring += f"    {dtype} r_pg = rate_eval.screened_rates(k_{self.secondary_rates[1].fname});\n"
            fstring += f"    {dtype} r_pa = rate_eval.screened_rates(k_{self.secondary_reverse[1].fname});\n"

            # now the approximation
            fstring += f"    {dtype} dd = 1.0_rt / (r_pg + r_pa);\n"
            fstring += "    rate = r_ag + r_ap * r_pg * dd;\n"
            fstring += "    if constexpr (std::is_same<T, rate_derivs_t>::value) {\n"
            fstring += f"        {dtype} drdT_ag = rate_eval.dscreened_rates_dT(k_{self.primary_rate.fname});\n"
            fstring += f"        {dtype} drdT_ap = rate_eval.dscreened_rates_dT(k_{self.secondary_rates[0].fname});\n"
            fstring += f"        {dtype} drdT_pg = rate_eval.dscreened_rates_dT(k_{self.secondary_rates[1].fname});\n"
            fstring += f"        {dtype} drdT_pa = rate_eval.dscreened_rates_dT(k_{self.secondary_reverse[1].fname});\n"
            fstring += "        drate_dT = drdT_ag + drdT_ap * r_pg * dd + r_ap * drdT_pg * dd - r_ap * r_pg * dd * dd * (drdT_pg + drdT_pa);\n"
            fstring += "    }\n"
        else:

            # first we need to get all of the rates that make this up
            fstring += f"    {dtype} r_ga = rate_eval.screened_rates(k_{self.primary_reverse.fname});\n"
            fstring += f"    {dtype} r_pa = rate_eval.screened_rates(k_{self.secondary_reverse[1].fname});\n"
            fstring += f"    {dtype} r_gp = rate_eval.screened_rates(k_{self.secondary_reverse[0].fname});\n"
            fstring += f"    {dtype} r_pg = rate_eval.screened_rates(k_{self.secondary_rates[1].fname});\n"

            # now the approximation
            fstring += f"    {dtype} dd = 1.0_rt / (r_pg + r_pa);\n"
            fstring += "    rate = r_ga + r_gp * r_pa * dd;\n"
            fstring += "    if constexpr (std::is_same<T, rate_derivs_t>::value) {\n"
            fstring += f"        {dtype} drdT_ga = rate_eval.dscreened_rates_dT(k_{self.primary_reverse.fname});\n"
            fstring += f"        {dtype} drdT_pa = rate_eval.dscreened_rates_dT(k_{self.secondary_reverse[1].fname});\n"
            fstring += f"        {dtype} drdT_gp = rate_eval.dscreened_rates_dT(k_{self.secondary_reverse[0].fname});\n"
            fstring += f"        {dtype} drdT_pg = rate_eval.dscreened_rates_dT(k_{self.secondary_rates[1].fname});\n"
            fstring += "        drate_dT = drdT_ga + drdT_gp * r_pa * dd + r_gp * drdT_pa * dd - r_gp * r_pa * dd * dd * (drdT_pg + drdT_pa);\n"
            fstring += "    }\n"

        if not leave_open:
            fstring += "}\n\n"

        return fstring
