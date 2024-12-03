"""
Classes and methods to interface with files storing rate data.
"""

import math
from pathlib import Path

import numpy as np

import pynucastro.numba_util as numba
from pynucastro.nucdata import Nucleus
from pynucastro.numba_util import jitclass
from pynucastro.rates.files import _find_rate_file


@jitclass([
    ('T9', numba.float64),
    ('T9i', numba.float64),
    ('T913', numba.float64),
    ('T913i', numba.float64),
    ('T953', numba.float64),
    ('lnT9', numba.float64)
])
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

    def __init__(self, T):
        """ return the Tfactors object.  Here, T is temperature in Kelvin """
        self.T9 = T/1.e9
        self.T9i = 1.0/self.T9
        self.T913i = self.T9i**(1./3.)
        self.T913 = self.T9**(1./3.)
        self.T953 = self.T9**(5./3.)
        self.lnT9 = np.log(self.T9)

    @property
    def array(self):
        """return t factors as array in order of lambda function"""
        return np.array([1, self.T9i, self.T913i, self.T913, self.T9, self.T953, self.lnT9])


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

        self.source = None
        self.modified = False

        # the fname is used when writing the code to evaluate the rate
        reactants_str = '_'.join([repr(nuc) for nuc in self.reactants])
        products_str = '_'.join([repr(nuc) for nuc in self.products])
        self.fname = f'{reactants_str}__{products_str}__{label}'

        if Q is None:
            self._set_q()
        else:
            self.Q = Q

        self.weak_type = weak_type

        # the identical particle factor scales the rate to prevent
        # double counting for a rate that has the same nucleus
        # multiple times as a reactant.  Usually we want this
        # behavior, but for approximate rates, sometimes we need to
        # disable it.
        self.use_identical_particle_factor = use_identical_particle_factor

        self._set_rhs_properties()
        self._set_screening()
        self._set_print_representation()

        self.tabular = False

        self.reverse = None

        self.rate_eval_needs_rho = False
        self.rate_eval_needs_comp = False

    def __repr__(self):
        return self.string

    def __hash__(self):
        return hash(self.__repr__())

    def __eq__(self, other):
        """ Determine whether two Rate objects are equal.
        They are equal if they contain identical reactants and products"""

        return (self.reactants, self.products) == (other.reactants, other.products)

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

    @property
    def id(self):
        return self.get_rate_id()

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

    def swap_protons(self):
        """Change any protons in the rate to NSE protons.  These
        act the same as protons but will be kept as distinct in
        the network."""

        p = Nucleus("p")
        p_nse = Nucleus("p_nse")

        for n, nuc in enumerate(self.reactants):
            if nuc == p:
                self.reactants[n] = p_nse

        for n, nuc in enumerate(self.products):
            if nuc == p:
                self.products[n] = p_nse

        # we need to update the Q value and the print string for the rate

        self._set_q()
        self._set_screening()
        self.fname = None    # reset so it will be updated
        self._set_print_representation()

    def modify_products(self, new_products):
        """
        change the products of the rate to new_products.  This will recompute
        the Q value and update the print respresentation.
        """

        self.products = Nucleus.cast_list(new_products, allow_single=True)
        self.modified = True

        # we need to update the Q value and the print string for the rate

        self._set_q()
        self._set_screening()
        self.fname = None    # reset so it will be updated
        self._set_print_representation()

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

    def eval(self, T, *, rho=None, comp=None):
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
            y_e_term = comp.ye
        else:
            y_e_term = 1.0

        # finally evaluate the rate -- for tabular rates, we need to set rhoY
        rate_eval = self.eval(T, rho=rho, comp=comp)

        return self.prefactor * dens_term * y_e_term * Y_term * rate_eval


class RateSource:
    """A class that stores the reference label information for various rates."""

    csv_path = _find_rate_file("rate_sources.csv")

    urls = {
        "debo": "https://doi.org/10.1103/RevModPhys.89.035007",
        "langanke": "https://doi.org/10.1006/adnd.2001.0865",
        "suzuki": "https://doi.org/10.3847/0004-637X/817/2/163",
        "reaclib": "https://reaclib.jinaweb.org/labels.php?action=viewLabel&label="
    }

    @staticmethod
    def _read_rate_sources(urls: dict[str, str], csv_path: Path) -> dict[str, dict[str, str]]:
        """Builds the labels dictionary from the supplied csv file."""

        labels = {}
        with csv_path.open("r") as csv:
            lines = csv.readlines()
            column_titles = lines[0].split("|")
            for line in lines[1:]:
                cells = [cell.strip() for cell in line.split("|")]
                label = cells[0]
                label_data = labels[label.lower()] = dict(zip(column_titles, cells))
                label_data["URL"] = urls.get(label, urls["reaclib"] + label)
        return labels

    labels = _read_rate_sources(urls, csv_path)

    @classmethod
    def source(cls, label: str) -> dict[str, str] | None:
        """Returns the source of a rate given its label, and None if not found."""

        return cls.labels.get(label.lower().strip())


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
