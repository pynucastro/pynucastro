"""
Classes and methods to interface with files storing rate data.
"""

import io
import math
import warnings
from collections import Counter
from enum import Enum
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

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


from pynucastro.constants import constants
from pynucastro.nucdata import Nucleus, UnsupportedNucleus

_pynucastro_dir = Path(__file__).parents[1]
_pynucastro_rates_dir = _pynucastro_dir/"library"
_pynucastro_tabular_dir = _pynucastro_rates_dir/"tabular"
_pynucastro_suzuki_dir = _pynucastro_tabular_dir/"suzuki"
_pynucastro_langanke_dir = _pynucastro_tabular_dir/"langanke"
_dirs = [
    _pynucastro_dir, _pynucastro_rates_dir, _pynucastro_tabular_dir,
    _pynucastro_suzuki_dir, _pynucastro_langanke_dir
]


def get_rates_dir() -> Path:
    return _pynucastro_rates_dir


def get_tabular_dir() -> Path:
    return _pynucastro_tabular_dir


class RateFileError(Exception):
    """An error occurred while trying to read a Rate from a file."""


def load_rate(rfile: str = None):
    """Try to load a rate of any type.

    :raises: :class:`.RateFileError`, :class:`.UnsupportedNucleus`
    """

    rate: Rate
    try:
        rate = TabularRate(rfile=rfile)
    except (RateFileError, UnsupportedNucleus):
        rate = ReacLibRate(rfile=rfile)

    return rate


def _find_rate_file(ratename: str | Path) -> Path:
    """locate the Reaclib or tabular rate or library file given its name.  Return
    None if the file cannot be located, otherwise return its path."""

    # check to see if the rate file is in the working dir,
    # is already the full path, or is in _dirs

    for path in ("", *_dirs):
        x = Path(path, ratename).resolve()
        if x.is_file():
            return x.resolve()

    # notify user we can't find the file
    raise RateFileError(f'File {ratename!r} not found in the working directory, {_pynucastro_rates_dir}, or {_pynucastro_tabular_dir}')


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
    :var T913i: 1.0 / T9 ** (1/3)
    :var T913:  T9 ** (1/3)
    :var T953:  T9 ** (5/3)
    :var lnT9:  log(T9)
    """

    def __init__(self, T: float) -> None:
        """ return the Tfactors object.  Here, T is temperature in Kelvin """
        self.T9 = T/1.e9
        self.T9i = 1.0/self.T9
        self.T913i = self.T9i**(1./3.)
        self.T913 = self.T9**(1./3.)
        self.T953 = self.T9**(5./3.)
        self.lnT9 = np.log(self.T9)

    @property
    def array(self) -> np.ndarray:
        """return t factors as array in order of lambda function"""
        return np.array([1, self.T9i, self.T913i, self.T913, self.T9, self.T953, self.lnT9])


class SingleSet:
    """ a set in Reaclib is one piece of a rate, in the form

        lambda = exp[ a_0 + sum_{i=1}^5  a_i T_9**(2i-5)/3  + a_6 log T_9]

    A single rate in Reaclib can be composed of multiple sets

    :param a: the coefficients of the exponential fit
    :param labelprops: a collection of flags that classify a ReacLib rate

    """

    def __init__(self, a, labelprops) -> None:
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

    def _update_label_properties(self) -> None:
        """ Set label and flags indicating Set is resonant,
            weak, or reverse. """
        assert isinstance(self.labelprops, str)
        assert len(self.labelprops) == 6

        self.label = self.labelprops[0:4]
        self.resonant = self.labelprops[4] == 'r'
        self.weak = self.labelprops[4] == 'w'
        self.reverse = self.labelprops[5] == 'v'

    def __eq__(self, other) -> bool:
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
        return lambda tf: float(np.exp(self.a[0] +
                                       self.a[1]*tf.T9i +
                                       self.a[2]*tf.T913i +
                                       self.a[3]*tf.T913 +
                                       self.a[4]*tf.T9 +
                                       self.a[5]*tf.T953 +
                                       self.a[6]*tf.lnT9))

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
        if all(q == 0.0 for q in self.a[1:]):
            string += "\namrex::ignore_unused(tfactors);"
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
    def __init__(self, reactants=None, products=None,
                 Q=None, weak_type="", label="generic",
                 use_identical_particle_factor=True):
        """a generic Rate class that acts as a base class for specific
        sources.  Here we only specify the reactants and products and Q value"""

        if reactants:
            self.reactants = Nucleus.cast_list(reactants)
        else:
            self.reactants = []

        if products:
            self.products = Nucleus.cast_list(products)
        else:
            self.products = []

        self.label = label

        # the fname is used when writing the code to evaluate the rate
        reactants_str = '_'.join([repr(nuc) for nuc in self.reactants])
        products_str = '_'.join([repr(nuc) for nuc in self.products])
        self.fname = f'{reactants_str}__{products_str}__{label}'

        if Q is None:
            self._set_q()
        else:
            self.Q = Q

        self.weak_type = weak_type

        self._set_rhs_properties()
        self._set_screening()
        self._set_print_representation()

        self.tabular = False

        self.reverse = None

        # the identical particle factor scales the rate to prevent
        # double counting for a rate that has the same nucleus
        # multiple times as a reactant.  Usually we want this
        # behavior, but for approximate rates, sometimes we need to
        # disable it.
        self.use_identical_particle_factor = use_identical_particle_factor

    def __repr__(self):
        return self.string

    def __hash__(self):
        return hash(self.__repr__())

    def __eq__(self, other):
        """ Determine whether two Rate objects are equal.
        They are equal if they contain identical reactants and products"""
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

        # from the masses of the nuclei, Q = M_products - M_reactants

        self.Q = 0
        for n in self.reactants:
            self.Q += n.mass
        for n in self.products:
            self.Q += -n.mass

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
        if self.use_identical_particle_factor:
            for r in set(self.reactants):
                self.inv_prefactor = self.inv_prefactor * math.factorial(self.reactants.count(r))
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

    def cname(self):
        """a C++-safe version of the rate name"""
        # replace the "__" separating reactants and products with "_to_"
        # and convert all other "__" to single "_"
        return self.fname.replace("__", "_to_", 1).replace("__", "_")

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

        ydot_string_components = []

        # prefactor
        if self.prefactor != 1.0:
            ydot_string_components.append(f"{self.prefactor:1.14e}")

        # density dependence
        if self.dens_exp == 1:
            ydot_string_components.append("rho")
        elif self.dens_exp != 0:
            ydot_string_components.append(f"rho**{self.dens_exp}")

        # electron fraction dependence
        if self.weak_type == 'electron_capture' and not self.tabular:
            ydot_string_components.append("ye(Y)")

        # composition dependence
        for r in sorted(set(self.reactants)):
            c = self.reactants.count(r)
            if c > 1:
                ydot_string_components.append(f"Y[j{r.raw}]**{c}")
            else:
                ydot_string_components.append(f"Y[j{r.raw}]")

        # rate_eval.{fname}
        ydot_string_components.append(f"rate_eval.{self.fname}")

        return "*".join(ydot_string_components)

    def eval(self, T, rhoY=None):
        raise NotImplementedError("base Rate class does not know how to eval()")

    def jacobian_string_py(self, y_i):
        """
        Return a string containing the term in a jacobian matrix
        in a reaction network corresponding to this rate differentiated
        with respect to y_i

        y_i is an object of the class ``Nucleus``.
        """
        if y_i not in self.reactants:
            return ""

        jac_string_components = []

        # prefactor
        if self.prefactor != 1.0:
            jac_string_components.append(f"{self.prefactor:1.14e}")

        # density dependence
        if self.dens_exp == 1:
            jac_string_components.append("rho")
        elif self.dens_exp != 0:
            jac_string_components.append(f"rho**{self.dens_exp}")

        # electron fraction dependence
        if self.weak_type == 'electron_capture' and not self.tabular:
            jac_string_components.append("ye(Y)")

        # composition dependence
        for r in sorted(set(self.reactants)):
            c = self.reactants.count(r)
            if y_i == r:
                # take the derivative
                if c == 1:
                    continue
                if c > 2:
                    jac_string_components.append(f"{c}*Y[j{r.raw}]**{c-1}")
                elif c == 2:
                    jac_string_components.append(f"2*Y[j{r.raw}]")
            else:
                # this nucleus is in the rate form, but we are not
                # differentiating with respect to it
                if c > 1:
                    jac_string_components.append(f"Y[j{r.raw}]**{c}")
                else:
                    jac_string_components.append(f"Y[j{r.raw}]")

        # rate_eval.{fname}
        jac_string_components.append(f"rate_eval.{self.fname}")

        return "*".join(jac_string_components)

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
        for r in sorted(set(self.reactants)):
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


class TableIndex(Enum):
    """a simple enum-like container for indexing the electron-capture tables"""
    RHOY = 0
    T = 1
    MU = 2
    DQ = 3
    VS = 4
    RATE = 5
    NU = 6
    GAMMA = 7


class ReacLibRate(Rate):
    """A single reaction rate.  Currently, this is a ReacLib rate, which
    can be composed of multiple sets, or a tabulated electron capture
    rate.

    :raises: :class:`.RateFileError`, :class:`.UnsupportedNucleus`
    """
    def __init__(self, rfile=None, chapter=None, original_source=None,
                 reactants=None, products=None, sets=None, labelprops=None, Q=None):
        """ rfile can be either a string specifying the path to a rate file or
        an io.StringIO object from which to read rate information. """
        # pylint: disable=super-init-not-called

        self.rfile_path = None
        self.rfile = None

        if isinstance(rfile, (str, Path)):
            rfile = Path(rfile)
            self.rfile_path = _find_rate_file(rfile)
            self.rfile = rfile.name

        self.chapter = chapter    # the Reaclib chapter for this reaction
        self.original_source = original_source   # the contents of the original rate file
        self.fname = None

        if reactants:
            self.reactants = Nucleus.cast_list(reactants)
        else:
            self.reactants = []

        if products:
            self.products = Nucleus.cast_list(products)
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

        self.use_identical_particle_factor = True

        if isinstance(rfile, Path):
            # read in the file, parse the different sets and store them as
            # SingleSet objects in sets[]
            f = self.rfile_path.open()
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
        self.products = Nucleus.cast_list(new_products, allow_single=True)
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

        x = (self.chapter == other.chapter) and (self.products == other.products) and \
                (self.reactants == other.reactants)
        if not x:
            return x
        x = len(self.sets) == len(other.sets)
        if not x:
            return x

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
                if check_chapter != self.chapter:
                    raise RateFileError(f'read chapter {check_chapter}, expected chapter {self.chapter} for this rate set.')
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

                chapter_dict = {
                    1: ((1,), (2,)),  # e1 -> e2
                    2: ((1,), (2, 3)),  # e1 -> e2 + e3
                    3: ((1,), (2, 3, 4)),  # e1 -> e2 + e3 + e4
                    4: ((1, 2), (3,)),  # e1 + e2 -> e3
                    5: ((1, 2), (3, 4)),  # e1 + e2 -> e3 + e4
                    6: ((1, 2), (3, 4, 5)),  # e1 + e2 -> e3 + e4 + e5
                    7: ((1, 2), (3, 4, 5, 6)),  # e1 + e2 -> e3 + e4 + e5 + e6
                    8: ((1, 2, 3), (4,)),  # e1 + e2 + e3 -> e4
                    9: ((1, 2, 3), (4, 5)),  # e1 + e2 + e3 -> e4 + e5
                    10: ((1, 2, 3, 4), (5, 6)),  # e1 + e2 + e3 + e4 -> e5 + e6
                    11: ((1,), (2, 3, 4, 5))  # e1 -> e2 + e3 + e4 + e5
                }

                try:
                    r, p = chapter_dict[self.chapter]
                    self.reactants += [Nucleus.from_cache(f[i-1]) for i in r]
                    self.products += [Nucleus.from_cache(f[j-1]) for j in p]

                    # support historical format, where chapter 8 also handles what are
                    # now chapter 9 rates
                    if self.chapter == 8 and len(f) == 5:
                        self.products.append(Nucleus.from_cache(f[4]))

                except KeyError as exc:
                    raise RateFileError(f'Chapter {self.chapter} could not be identified in {self.original_source}') from exc

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

    def write_to_file(self, f) -> None:
        """ Given a file object, write rate data to the file. """

        if self.original_source is None:
            raise NotImplementedError(
                f"Original source is not stored for this rate ({self})."
                " At present, we cannot reconstruct the rate representation without"
                " storing the original source."
            )

        print(self.original_source, file=f)

    def get_rate_id(self) -> str:
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

    def function_string_cxx(self, dtype="double", specifiers="inline", leave_open=False, extra_args=()):
        """
        Return a string containing C++ function that computes the
        rate
        """

        args = ["const tf_t& tfactors", f"{dtype}& rate", f"{dtype}& drate_dT", *extra_args]
        fstring = ""
        fstring += "template <int do_T_derivatives>\n"
        fstring += f"{specifiers}\n"
        fstring += f"void rate_{self.cname()}({', '.join(args)}) {{\n\n"
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
        _ = rhoY  # unused by this subclass

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
        _ = (rhoYmin, rhoYmax)  # unused by this subclass

        fig, ax = plt.subplots(figsize=figsize)

        temps = np.logspace(np.log10(Tmin), np.log10(Tmax), 100)
        r = np.zeros_like(temps)

        for n, T in enumerate(temps):
            r[n] = self.eval(T)

        ax.loglog(temps, r)
        ax.set_xlabel(r"$T$")

        if self.dens_exp == 0:
            ax.set_ylabel(r"$\tau$")
        elif self.dens_exp == 1:
            ax.set_ylabel(r"$N_A <\sigma v>$")
        elif self.dens_exp == 2:
            ax.set_ylabel(r"$N_A^2 <n_a n_b n_c v>$")

        ax.set_title(fr"{self.pretty_string}")

        return fig


if numba is not None:
    interpolator_spec = [
        ('data', numba.float64[:, :]),
        ('table_rhoy_lines', numba.int32),
        ('table_temp_lines', numba.int32),
        ('rhoy', numba.float64[:]),
        ('temp', numba.float64[:])
    ]
else:
    interpolator_spec = []


@jitclass(interpolator_spec)
class TableInterpolator:
    """A simple class that holds a pointer to the table data and
    methods that allow us to interpolate a variable"""

    def __init__(self, table_rhoy_lines, table_temp_lines, table_data):

        self.data = table_data
        self.table_rhoy_lines = table_rhoy_lines
        self.table_temp_lines = table_temp_lines

        # for easy indexing, store a 1-d array of T and rhoy
        self.rhoy = self.data[::self.table_temp_lines, TableIndex.RHOY.value]
        self.temp = self.data[0:self.table_temp_lines, TableIndex.T.value]

    def _get_logT_idx(self, logt0):
        """return the index into the temperatures such that
        T[i-1] < t0 <= T[i].  We return i-1 here, corresponding to
        the lower value.
        Note: we work in terms of log10()
        """

        max_idx = len(self.temp) - 1
        return max(0, min(max_idx, np.searchsorted(self.temp, logt0)) - 1)

    def _get_logrhoy_idx(self, logrhoy0):
        """return the index into rho*Y such that
        rhoY[i-1] < rhoy0 <= rhoY[i].  We return i-1 here,
        corresponding to the lower value.
        Note: we work in terms of log10()

        """

        max_idx = len(self.rhoy) - 1
        return max(0, min(max_idx, np.searchsorted(self.rhoy, logrhoy0)) - 1)

    def _rhoy_T_to_idx(self, irhoy, jtemp):
        """given a pair (irhoy, jtemp) into the table, return the 1-d index
        into the underlying data array assuming row-major ordering"""

        return irhoy * self.table_temp_lines + jtemp

    def interpolate(self, logrhoy, logT, component):
        """given logrhoy and logT, do bilinear interpolation to
        find the value of the data component in the table"""

        # We are going to do bilinear interpolation.  We create a
        # polynomial of the form:
        #
        # f = A [log(rho) - log(rho_i)] [log(T) - log(T_j)] +
        #     B [log(rho) - log(rho_i)] +
        #     C [log(T) - log(T_j)] +
        #     D
        #
        # we then find the i,j such that our point is in the
        # box with corners (i,j) to (i+1,j+1), and solve for
        # A, B, C, D

        # find the T and rhoY in the data table corresponding to the
        # lower left

        if logT < self.temp.min() or logT > self.temp.max():
            raise ValueError("temperature out of table bounds")

        if logrhoy < self.rhoy.min() or logrhoy > self.rhoy.max():
            raise ValueError("rhoy out of table bounds")

        irhoy = self._get_logrhoy_idx(logrhoy)
        jT = self._get_logT_idx(logT)

        # note: rhoy and T are already stored as log

        dlogrho = self.rhoy[irhoy+1] - self.rhoy[irhoy]
        dlogT = self.temp[jT+1] - self.temp[jT]

        # get the data at the 4 points

        idx = self._rhoy_T_to_idx(irhoy, jT)
        f_ij = self.data[idx, component]

        idx = self._rhoy_T_to_idx(irhoy+1, jT)
        f_ip1j = self.data[idx, component]

        idx = self._rhoy_T_to_idx(irhoy, jT+1)
        f_ijp1 = self.data[idx, component]

        idx = self._rhoy_T_to_idx(irhoy+1, jT+1)
        f_ip1jp1 = self.data[idx, component]

        D = f_ij
        C = (f_ijp1 - f_ij) / dlogT
        B = (f_ip1j - f_ij) / dlogrho
        A = (f_ip1jp1 - B * dlogrho - C * dlogT - D) / (dlogrho * dlogT)

        r = (A * (logrhoy - self.rhoy[irhoy]) * (logT - self.temp[jT]) +
             B * (logrhoy - self.rhoy[irhoy]) + C * (logT - self.temp[jT]) + D)

        return r


class TabularRate(Rate):
    """A tabular rate.

    :raises: :class:`.RateFileError`, :class:`.UnsupportedNucleus`
    """
    def __init__(self, rfile=None):
        """ rfile can be either a string specifying the path to a rate file or
        an io.StringIO object from which to read rate information. """
        super().__init__()

        self.rfile_path = None
        self.rfile = None

        if isinstance(rfile, (str, Path)):
            rfile = Path(rfile)
            self.rfile_path = _find_rate_file(rfile)
            self.rfile = rfile.name

        self.fname = None

        self.label = "tabular"
        self.tabular = True

        # we should initialize this somehow
        self.weak_type = ""

        if isinstance(rfile, Path):
            # read in the file, parse the different sets and store them as
            # SingleSet objects in sets[]
            f = self.rfile_path.open()
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

        # store the extrema of the thermodynamics
        _rhoy = self.tabular_data_table[::self.table_temp_lines, TableIndex.RHOY.value]
        _temp = self.tabular_data_table[0:self.table_temp_lines, TableIndex.T.value]

        self.table_Tmin = 10.0**(_temp.min())
        self.table_Tmax = 10.0**(_temp.max())
        self.table_rhoYmin = 10.0**(_rhoy.min())
        self.table_rhoYmax = 10.0**(_rhoy.max())

        self.interpolator = TableInterpolator(self.table_rhoy_lines, self.table_temp_lines,
                                              self.tabular_data_table)

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
        if self.chapter != "t":
            raise RateFileError(f"Invalid chapter for TabularRate ({self.chapter})")

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
        except UnsupportedNucleus as ex:
            raise RateFileError(f'Nucleus objects could not be identified in {self.original_source}') from ex

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

        # since the reactants and products were only now set, we need
        # to recompute Q -- this is used for finding rate pairs
        self._set_q()

    def _set_rhs_properties(self):
        """ compute statistical prefactor and density exponent from the reactants. """
        self.prefactor = 1.0  # this is 1/2 for rates like a + a (double counting)
        self.inv_prefactor = 1
        if self.use_identical_particle_factor:
            for r in set(self.reactants):
                self.inv_prefactor = self.inv_prefactor * math.factorial(self.reactants.count(r))
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

        fstring += f"    {self.fname}_interpolator = TableInterpolator(*{self.fname}_info)\n"

        fstring += f"    r = {self.fname}_interpolator.interpolate(np.log10(rhoY), np.log10(T), TableIndex.RATE.value)\n"
        fstring += f"    rate_eval.{self.fname} = 10.0**r\n\n"

        return fstring

    def get_tabular_rate(self):
        """read the rate data from .dat file """

        # find .dat file and read it
        self.table_path = _find_rate_file(self.table_file)
        t_data2d = []
        with self.table_path.open() as tabular_file:
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
        self.tabular_data_table = np.array(t_data2d, dtype=np.float64)

    def eval(self, T, rhoY=None):
        """ evauate the reaction rate for temperature T """

        r = self.interpolator.interpolate(np.log10(rhoY), np.log10(T),
                                          TableIndex.RATE.value)
        return 10.0**r

    def get_nu_loss(self, T, rhoY):
        """ get the neutrino loss rate for the reaction if tabulated"""

        r = self.interpolator.interpolate(np.log10(rhoY), np.log10(T),
                                          TableIndex.NU.value)
        return 10**r

    def plot(self, *, Tmin=None, Tmax=None, rhoYmin=None, rhoYmax=None,
             color_field='rate', figsize=(10, 10)):
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

        if Tmin is None:
            Tmin = self.table_Tmin
        if Tmax is None:
            Tmax = self.table_Tmax
        if rhoYmin is None:
            rhoYmin = self.table_rhoYmin
        if rhoYmax is None:
            rhoYmax = self.table_rhoYmax

        data = self.tabular_data_table

        inde1 = data[:, TableIndex.T.value] <= np.log10(Tmax)
        inde2 = data[:, TableIndex.T.value] >= np.log10(Tmin)
        inde3 = data[:, TableIndex.RHOY.value] <= np.log10(rhoYmax)
        inde4 = data[:, TableIndex.RHOY.value] >= np.log10(rhoYmin)
        data_heatmap = data[inde1 & inde2 & inde3 & inde4].copy()

        rows, row_pos = np.unique(data_heatmap[:, 0], return_inverse=True)
        cols, col_pos = np.unique(data_heatmap[:, 1], return_inverse=True)
        pivot_table = np.zeros((len(rows), len(cols)), dtype=data_heatmap.dtype)

        if color_field == 'rate':
            icol = TableIndex.RATE.value
            title = f"{self.weak_type} rate in log10(1/s)"
            cmap = 'magma'

        elif color_field == 'nu_loss':
            icol = TableIndex.NU.value
            title = "neutrino energy loss rate in log10(erg/s)"
            cmap = 'viridis'

        else:
            raise ValueError("color_field must be either 'rate' or 'nu_loss'.")

        try:
            pivot_table[row_pos, col_pos] = data_heatmap[:, icol]
        except ValueError:
            print("Divide by zero encountered in log10\nChange the scale of T or rhoY")

        im = ax.imshow(pivot_table, cmap=cmap, origin="lower",
                       extent=[np.log10(Tmin), np.log10(Tmax), np.log10(rhoYmin), np.log10(rhoYmax)])
        fig.colorbar(im, ax=ax)

        ax.set_xlabel(r"$\log(T)$ [K]")
        ax.set_ylabel(r"$\log(\rho Y_e)$ [g/cm$^3$]")
        ax.set_title(fr"{self.pretty_string}" + "\n" + title)

        return fig


class DerivedRate(ReacLibRate):
    """
    This class is a derived class from `Rate` with the purpose of computing the inverse rate
    by the application of detailed balance to the forward reactions.
    """

    def __init__(self, rate, compute_Q=False, use_pf=False):

        self.use_pf = use_pf
        self.rate = rate
        self.compute_Q = compute_Q

        if not isinstance(rate, Rate):
            raise TypeError('rate must be a Rate subclass')

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
            prefactor += -np.log(constants.N_A) * (len(self.rate.reactants) - len(self.rate.products))

            for nucr in self.rate.reactants:
                prefactor += 1.5*np.log(nucr.A) + np.log(nucr.spin_states)
                Q += nucr.A_nuc
            for nucp in self.rate.products:
                prefactor += -1.5*np.log(nucp.A) - np.log(nucp.spin_states)
                Q -= nucp.A_nuc

            if self.compute_Q:
                Q = Q * constants.m_u_MeV
            else:
                Q = self.rate.Q

            prefactor += np.log(self.counter_factors()[1]) - np.log(self.counter_factors()[0])

            if len(self.rate.reactants) == len(self.rate.products):
                prefactor += 0.0
            else:
                F = (constants.m_u * constants.k * 1.0e9 / (2.0*np.pi*constants.hbar**2))**(1.5*(len(self.rate.reactants) - len(self.rate.products)))
                prefactor += np.log(F)

            a_rev = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            a_rev[0] = prefactor + a[0]
            a_rev[1] = a[1] - Q / (1.0e9 * constants.k_MeV)
            a_rev[2] = a[2]
            a_rev[3] = a[3]
            a_rev[4] = a[4]
            a_rev[5] = a[5]
            a_rev[6] = a[6] + 1.5*(len(self.rate.reactants) - len(self.rate.products))
            sset_d = SingleSet(a=a_rev, labelprops=rate.labelprops)
            derived_sets.append(sset_d)

        super().__init__(rfile=self.rate.rfile, chapter=self.rate.chapter, original_source=self.rate.original_source,
                reactants=self.rate.products, products=self.rate.reactants, sets=derived_sets, labelprops="derived", Q=-Q)

    def _warn_about_missing_pf_tables(self):
        skip_nuclei = {Nucleus("h1"), Nucleus("n"), Nucleus("he4")}
        for nuc in set(self.rate.reactants + self.rate.products) - skip_nuclei:
            if not nuc.partition_function:
                warnings.warn(UserWarning(f'{nuc} partition function is not supported by tables: set pf = 1.0 by default'))

    def eval(self, T, rhoY=None):

        r = super().eval(T=T, rhoY=rhoY)
        z_r = 1.0
        z_p = 1.0
        if self.use_pf:
            self._warn_about_missing_pf_tables()

            for nucr in self.rate.reactants:
                if not nucr.partition_function:
                    continue
                    #nucr.partition_function = lambda T: 1.0
                z_r *= nucr.partition_function.eval(T)

            for nucp in self.rate.products:
                if not nucp.partition_function:
                    continue
                    #nucp.partition_function = lambda T: 1.0
                z_p *= nucp.partition_function.eval(T)

            return r*z_r/z_p
        return r

    def function_string_py(self):
        """
        Return a string containing python function that computes the
        rate
        """

        self._warn_about_missing_pf_tables()

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

    def function_string_cxx(self, dtype="double", specifiers="inline", leave_open=False, extra_args=()):
        """
        Return a string containing C++ function that computes the
        rate
        """

        self._warn_about_missing_pf_tables()

        extra_args = ["[[maybe_unused]] part_fun::pf_cache_t& pf_cache", *extra_args]
        fstring = super().function_string_cxx(dtype=dtype, specifiers=specifiers, leave_open=True, extra_args=extra_args)

        # right now we have rate and drate_dT without the partition function
        # now the partition function corrections

        if self.use_pf:

            fstring += "\n"
            for nuc in set(self.rate.reactants + self.rate.products):
                fstring += f"    {dtype} {nuc}_pf, d{nuc}_pf_dT;\n"

                if nuc.partition_function:
                    fstring += f"    // interpolating {nuc} partition function\n"
                    fstring += f"    get_partition_function_cached({nuc.cindex()}, tfactors, pf_cache, {nuc}_pf, d{nuc}_pf_dT);\n"
                else:
                    fstring += f"    // setting {nuc} partition function to 1.0 by default, independent of T\n"
                    fstring += f"    {nuc}_pf = 1.0_rt;\n"
                    fstring += f"    d{nuc}_pf_dT = 0.0_rt;\n"
                fstring += "\n"

            fstring += f"    {dtype} z_r = "
            fstring += " * ".join([f"{nucr}_pf" for nucr in self.rate.reactants])
            fstring += ";\n"

            fstring += f"    {dtype} z_p = "
            fstring += " * ".join([f"{nucp}_pf" for nucp in self.rate.products])
            fstring += ";\n\n"

            # now the derivatives, via chain rule
            chain_terms = []
            for n in self.rate.reactants:
                chain_terms.append(" * ".join([f"{nucr}_pf" for nucr in self.rate.reactants if nucr != n] + [f"d{n}_pf_dT"]))

            fstring += f"    {dtype} dz_r_dT = "
            fstring += " + ".join(chain_terms)
            fstring += ";\n"

            chain_terms = []
            for n in self.rate.products:
                chain_terms.append(" * ".join([f"{nucp}_pf" for nucp in self.rate.products if nucp != n] + [f"d{n}_pf_dT"]))

            fstring += f"    {dtype} dz_p_dT = "
            fstring += " + ".join(chain_terms)
            fstring += ";\n\n"

            fstring += f"    {dtype} dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);\n\n"

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
            reactant_factor *= math.factorial(react_counts[nuc])

        product_factor = 1.0
        for nuc in set(self.rate.products):
            product_factor *= math.factorial(prod_counts[nuc])

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
