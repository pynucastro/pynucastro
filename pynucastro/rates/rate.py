"""Classes and methods to interface with files storing rate data."""

import math
from pathlib import Path

import numpy as np

import pynucastro.numba_util as numba
from pynucastro.nucdata import Nucleus
from pynucastro.numba_util import jitclass
from pynucastro.rates.files import _find_rate_file
from pynucastro.screening import (get_screening_map, make_plasma_state,
                                  make_screen_factors)


class BaryonConservationError(Exception):
    """Exception for the case where we don't have the same number of
    baryons on the left and righthand sides of the reaction.

    """


@jitclass([
    ('T9', numba.float64),
    ('T9i', numba.float64),
    ('T913', numba.float64),
    ('T913i', numba.float64),
    ('T953', numba.float64),
    ('lnT9', numba.float64)
])
class Tfactors:
    """Precompute temperature factors for speed

    Parameters
    ----------
    T : float
        Temperature in Kelvin.

    Attributes
    ----------
    T9 : float
        Temperature divided by 1.e9 K.
    T9i : float
        1.0 / T9
    T913i : float
        1.0 / T9**(1/3)
    T913 : float
        T9**(1/3)
    T953 : float
        T9**(5/3)
    lnT9 : float
        log(T9)

    """

    def __init__(self, T):
        self.T9 = T/1.e9
        self.T9i = 1.0/self.T9
        self.T913i = self.T9i**(1./3.)
        self.T913 = self.T9**(1./3.)
        self.T953 = self.T9**(5./3.)
        self.lnT9 = np.log(self.T9)

    @property
    def array(self):
        """Compute the various T factors in the same order required by
        the lambda fit used by ReacLib rates

        Returns
        -------
        numpy.ndarray

        """
        return np.array([1, self.T9i, self.T913i, self.T913,
                         self.T9, self.T953, self.lnT9])


class Rate:
    """The base reaction rate class.  Most rate types will subclass
    this and extend to their particular format.

    Parameters
    ----------
    reactants : list(str), list(Nucleus)
        the reactants for the reaction
    products : list(str), list(Nucleus)
        the products of the reactions
    Q : float
        the energy release (in MeV) for the rate
    weak_type : str
        the type of weak reaction the rate represents.
        Possible values include "electron_capture", "beta_decay", or
        may include "_pos_" or "_neg_".
    label : str
        A descriptive label for the rate (usually representative of the
        source
    use_identical_particle_factor : bool
        For a rate that has the same nucleus multiple times as a reactant
        we apply a multiplicity factor, N!, where N is the number of times
        the nucleus appears.  This can be disabled by setting
        use_identical_particle_factor = False

    """

    def __init__(self, reactants=None, products=None,
                 Q=None, weak_type="", label="generic",
                 stoichiometry=None, rate_source=None,
                 use_identical_particle_factor=True):

        if reactants:
            self.reactants = Nucleus.cast_list(reactants)
        else:
            self.reactants = []

        if products:
            self.products = Nucleus.cast_list(products)
        else:
            self.products = []

        self.label = label

        if rate_source is not None:
            self.source = RateSource.source(rate_source)

        self.weak = False
        self.weak_type = weak_type
        if self.weak_type:
            self.weak = True

        self.tabular = False
        self.derived_from_inverse = False

        # the identical particle factor scales the rate to prevent
        # double counting for a rate that has the same nucleus
        # multiple times as a reactant.  Usually we want this
        # behavior, but for approximate rates, sometimes we need to
        # disable it.
        self.use_identical_particle_factor = use_identical_particle_factor

        # some subclasses might define a stoichmetry as a dict{Nucleus}
        # that gives the numbers for the dY/dt equations
        self.stoichiometry = stoichiometry

        # Set Q-value of the reaction rate. Needs to go after stoichiometry.
        if Q is None:
            self._set_q()
        else:
            self.Q = Q

        self._set_rhs_properties()
        self._set_screening()
        self._set_print_representation()

        # ensure that baryon number is conserved
        test = (
            sum(n.A * self.reactant_count(n) for n in set(self.reactants)) ==
            sum(n.A * self.product_count(n) for n in set(self.products))
        )

        if not test:
            raise BaryonConservationError(f"baryon number not conserved in rate {self}")

        self.rate_eval_needs_rho = False
        self.rate_eval_needs_comp = False

    def __repr__(self):
        return self.string

    def __hash__(self):
        return hash(self.__repr__())

    def __eq__(self, other):
        """Determine whether two Rate objects are equal.  They are
        equal if they contain identical reactants and products

        """

        return (self.reactants, self.products) == (other.reactants, other.products)

    def __lt__(self, other):
        """Sort such that lightest reactants come first, and then look
        at products.

        """

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
        """Set the Q value of the reaction (in MeV)."""

        # from the masses of the nuclei, Q = M_products - M_reactants

        self.Q = 0
        for n in set(self.reactants):
            c = self.reactant_count(n)
            self.Q += c * n.mass
        for n in set(self.products):
            c = self.product_count(n)
            self.Q += -c * n.mass

    def _set_print_representation(self):
        """Compose the string representations of this Rate.
        This includes string,rid, pretty_string, and fname.
        String is output to the terminal, rid is used as a dict key,
        and pretty_string is latex, and fname is used when writing the
        code to evaluate the rate.

        """

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

        reactant_Zs = sum(n.Z * self.reactant_count(n) for n in set(self.reactants))
        product_Zs = sum(n.Z * self.product_count(n) for n in set(self.products))

        reactant_As = sum(n.A * self.reactant_count(n) for n in set(self.reactants))
        product_As = sum(n.A * self.product_count(n) for n in set(self.products))

        strong_test = reactant_Zs == product_Zs and reactant_As == product_As

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

            elif self.weak_type and "_pos_" in self.weak_type:

                # we expect a positron on the right -- let's make sure
                assert sum(n.Z for n in self.reactants) == sum(n.Z for n in self.products) + 1

                rhs_other.append("e+")
                rhs_other.append("nu")

            elif self.weak_type and "_neg_" in self.weak_type:

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

        # this produces a sorted list with no dupes
        react_set = list(dict.fromkeys(treactants))
        for n, r in enumerate(react_set):
            c = self.reactant_count(r)
            if c == 2:
                # special case so we do C12 + C12 instead of 2 C12
                self.string += f"{r.c()} + {r.c()}"
                self.rid += f"{r} + {r}"
                self.pretty_string += fr"{r.pretty} + {r.pretty}"
            else:
                factor = ""
                if c != 1:
                    factor = f"{c} "
                self.string += f"{factor}{r.c()}"
                self.rid += f"{factor}{r}"
                if c != 1:
                    self.pretty_string += fr"{factor}~{r.pretty}"
                else:
                    self.pretty_string += fr"{r.pretty}"
            if not n == len(react_set)-1:
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

        prod_set = list(dict.fromkeys(self.products))
        for n, p in enumerate(prod_set):
            c = self.product_count(p)
            if c == 2:
                # special case for 2 species
                self.string += f"{p.c()} + {p.c()}"
                self.rid += f"{p} + {p}"
                self.pretty_string += fr"{p.pretty} + {p.pretty}"
            else:
                factor = ""
                if c != 1:
                    factor = f"{c} "
                self.string += f"{factor}{p.c()}"
                self.rid += f"{factor}{p}"
                self.pretty_string += fr"{factor}{p.pretty}"
            if not n == len(prod_set)-1:
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

        reactants_str = '_'.join([repr(nuc) for nuc in self.reactants])
        products_str = '_'.join([repr(nuc) for nuc in self.products])
        self.fname = f'{reactants_str}_to_{products_str}_{self.label}'

    def _set_rhs_properties(self):
        """Compute statistical prefactor and density exponent from the
        reactants.

        """
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
        """Determine if this rate is eligible for screening and the
        nuclei to use.

        """
        # Tells if this rate is eligible for screening, and if it is
        # then Rate.ion_screen is a 2-element (3 for 3-alpha) list of
        # Nucleus objects for screening; otherwise it is set to none
        self.ion_screen = []
        nucz = [q for q in self.reactants if q.Z != 0]
        if len(nucz) > 1:
            nucz.sort(key=lambda x: x.Z)
            self.ion_screen = []
            self.ion_screen.append(nucz[0])
            self.ion_screen.append(nucz[1])
            if len(nucz) == 3:
                self.ion_screen.append(nucz[2])

    def get_rate_id(self):
        """Get an identifying string for this rate.

        Returns
        -------
        str

        """
        return f'{self.rid} <{self.label.strip()}>'

    @property
    def id(self):
        """Get the rate's id string

        Returns
        -------
        str

        """
        return self.get_rate_id()

    def heaviest(self):
        """Get the heaviest nuclide in this Rate.  If two nuclei are
        tied in mass number, return the one with the lowest atomic
        number.

        Returns
        -------
        Nucleus

        """
        nuc = self.reactants[0]
        for n in self.reactants + self.products:
            if n.A > nuc.A or (n.A == nuc.A and n.Z < nuc.Z):
                nuc = n
        return nuc

    def lightest(self):
        """Get the lightest nuclide in this Rate.  If two nuclei are
        tied in mass number, return the one with the highest atomic
        number.

        Returns
        -------
        Nucleus

        """
        nuc = self.reactants[0]
        for n in self.reactants + self.products:
            if n.A < nuc.A or (n.A == nuc.A and n.Z > nuc.Z):
                nuc = n
        return nuc

    def swap_protons(self):
        """Change any protons in the rate to NSE protons.  These act
        the same as protons but will be kept as distinct in the
        network.

        """

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
        self._set_print_representation()

    def reactant_count(self, n):
        """Return the number of times nucleus n appears as a reactant
        in the rate.  Use the stoichiometry dict if present.

        Parameters
        ----------
        n : Nucleus
            the nucleus appearing as a reactant

        Returns
        -------
        float

        """

        c_reac = self.reactants.count(n)

        if self.stoichiometry and c_reac > 0:
            return self.stoichiometry.get(n, c_reac)
        return c_reac

    def product_count(self, n):
        """Return the number of times nucleus n appears as a product
        in the rate.  Use the stoichiometry dict if present.

        Parameters
        ----------
        n : Nucleus
            the nucleus appearing as a product

        Returns
        -------
        float

        """

        c_prod = self.products.count(n)

        if self.stoichiometry and c_prod > 0:
            return self.stoichiometry.get(n, c_prod)
        return c_prod

    def modify_products(self, new_products):
        """Change the products of the rate to new_products.  This will
        recompute the Q value and update the print representation.

        Parameters
        ----------
        new_products : list(Nucleus)
            the new products to use with the rate.

        """

        self.products = Nucleus.cast_list(new_products, allow_single=True)
        self.modified = True

        # we need to update the Q value and the print string for the rate

        self._set_q()
        self._set_screening()
        self._set_print_representation()

    def evaluate_screening(self, rho, T, composition, screen_func):
        """Evaluate the screening correction for this rate.

        Parameters
        ----------
        rho : float
            density used to evaluate screening
        T : float
            temperature used to evaluate screening
        composition : Composition
            composition used to evaluate screening
        screen_func : Callable
            one of the screening functions from :py:mod:`pynucastro.screening`

        Returns
        -------
        float

        """

        ys = composition.get_molar()
        plasma_state = make_plasma_state(T, rho, ys)
        scor = 1.0

        # We can have three cases:
        # 3-body reaction, i.e. 3-alpha: 2 ScreeningPair's
        # 2-body reaction              : 1 ScreeningPair
        # Photodisintegration (1-body) : 0 ScreeningPair

        screening_map = get_screening_map([self])

        # Handle 0 ScreeningPair case
        if not screening_map:
            return scor

        # Handle 3-alpha case explicitly
        if "He4_He4_He4" in [scr.name for scr in screening_map]:

            # We should have two ScreeningPair's in this case
            # He4_He4_He4 and He4_He4_He4_dummy
            assert len(screening_map) == 2

            for scr in screening_map:
                scn_fac = make_screen_factors(scr.n1, scr.n2)
                scor *= screen_func(plasma_state, scn_fac)

        # Now handle 2-body reaction
        else:
            scr = screening_map[0]

            # Make sure no dummy nuclei exist in this case
            # Otherwise we have more than 2-body reaction
            assert not (scr.n1.dummy or scr.n2.dummy)
            assert len(screening_map) == 1

            scn_fac = make_screen_factors(scr.n1, scr.n2)
            scor = screen_func(plasma_state, scn_fac)

        return scor

    def ydot_string_py(self):
        """Construct the string containing the term in a dY/dt
        equation in a reaction network corresponding to this rate.

        Returns
        -------
        str

        """

        ydot_string_components = []

        # prefactor (for double counting)
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

    def eval(self, T, *, rho=None, comp=None,
             screen_func=None):
        """Evaluate the reaction rate for temperature T.  This is a stub
        and should be implemented by the derived class.

        Parameters
        ----------
        T : float
            the temperature to evaluate the rate at
        rho : float
            the density to evaluate the rate and screening effects at.
        comp : float
            the composition (of type
            :py:class:`Composition <pynucastro.networks.rate_collection.Composition>`)
            to evaluate the rate and screening effects with.
        screen_func : Callable
            one of the screening functions from :py:mod:`pynucastro.screening`
            -- if provided, then the rate will include the screening correction

        Raises
        ------
        NotImplementedError

        """

        raise NotImplementedError("base Rate class does not know how to eval()")

    def function_string_py(self):
        """Return a string containing the python function that
        computes the rate.

        Raises
        ------
        NotImplementedError

        """

        raise NotImplementedError("base Rate class does not implement function_string_py()")

    def jacobian_string_py(self, y_i):
        """Return a string containing the term in a jacobian matrix
        in a reaction network corresponding to this rate
        differentiated with respect to ``y_i``

        Parameters
        ----------
        y_i : Nucleus
            the nucleus we are differentiating with respect to.

        Returns
        -------
        str

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

    def eval_jacobian_term(self, T, rho, comp, y_i, *,
                           screen_func=None):
        """Evaluate drate/d(y_i), the derivative of the rate with
        respect to ``y_i``.  This rate term has the full composition
        and density dependence, i.e.:

        rate = ρ**n Y1**a Y2**b ... N_A <σv>

        The derivative is only non-zero if this term depends on
        nucleus ``y_i``.

        Parameters
        ----------
        T : float
            the temperature to evaluate the rate with
        rho : float
            the density to evaluate the rate with
        comp : Composition
            the composition to use in the rate evaluation
        y_i : Nucleus
            the nucleus we are differentiating with respect to
        screen_func : Callable
            one of the screening functions from :py:mod:`pynucastro.screening`
            -- if provided, then the jacobian_term will include the
            screening correction.

        Returns
        -------
        float

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
        rate_eval = self.eval(T, rho=rho, comp=comp,
                              screen_func=screen_func)

        return self.prefactor * dens_term * y_e_term * Y_term * rate_eval


class RateSource:
    """A class that stores the reference label information for various rates."""

    csv_path = _find_rate_file("rate_sources.csv")

    urls = {
        "debo": "https://doi.org/10.1103/RevModPhys.89.035007",
        "langanke": "https://doi.org/10.1006/adnd.2001.0865",
        "suzuki": "https://doi.org/10.3847/0004-637X/817/2/163",
        "ffn": "https://doi.org/10.1086/190779",
        "pruet_fuller": "https://doi.org/10.1086/376753",
        "reaclib": "https://reaclib.jinaweb.org/labels.php?action=viewLabel&label=",
        "iliadis2022": "https://journals.aps.org/prc/abstract/10.1103/PhysRevC.106.055802"
    }

    @staticmethod
    def _read_rate_sources(urls: dict[str, str], csv_path: Path) -> dict[str, dict[str, str]]:
        """Build the labels dictionary from the supplied csv file."""

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
        """Return the source of a rate given its label, and None if
        not found.

        Parameters
        ----------
        label : str

        """

        return cls.labels.get(label.lower().strip())


class RatePair:
    """A pair of rates: the forward and reverse rates for a single
    reaction sequence.  Forward rates are those with Q >= 0.

    Parameters
    ----------
    forward : Rate
        the forward rate
    reverse : Rate
        the reverse rate

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
