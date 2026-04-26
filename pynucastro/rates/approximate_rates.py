"""Classes and methods for describing approximations to rates.  Often
this involves using a group of existing rates to enforce an
equilibrium through a nucleus.

"""

import math

from pynucastro.nucdata import Nucleus
from pynucastro.rates.rate import Rate


def _assert_rate_prop(rate, *,
                      reactants=None, products=None,
                      num_reactants=None, num_products=None):

    if reactants:
        for nuc in reactants:
            assert nuc in rate.reactants
    if num_reactants:
        assert len(rate.reactants) == num_reactants

    if products:
        for nuc in products:
            assert nuc in rate.products
    if num_products:
        assert len(rate.products) == num_products


def create_double_neutron_capture(lib, reactant, product):
    """Return a pair of :py:class:`ApproximateRate` objects for the
    A(n,g)X(n,g)B -> A(nn,g)B approximation

    Parameters
    ----------
    lib : Library
         A Library object containing the neutron-capture rates
    reactant : Nucleus, str
         The reactant, A, in the sequence A(n,g)X(n,g)B
    product: Nucleus, str
         The product, B, in the sequence A(n,g)X(n,g)B

    Returns
    -------
    ApproximateRate, ApproximateRate

    """

    if isinstance(reactant, str):
        reactant = Nucleus(reactant)

    if isinstance(product, str):
        product = Nucleus(product)

    intermediate = reactant + Nucleus("n")

    rates = {}
    rates["A(n,g)X"] = lib.get_rate_by_name(f"{reactant.raw}(n,){intermediate.raw}")
    rates["X(n,g)B"] = lib.get_rate_by_name(f"{intermediate.raw}(n,){product.raw}")

    rates["B(g,n)X"] = lib.get_rate_by_name(f"{product.raw}(,n){intermediate.raw}")
    rates["X(g,n)A"] = lib.get_rate_by_name(f"{intermediate.raw}(,n){reactant.raw}")

    forward = ApproximateRate(rates, approx_type="nn_g",
                              use_identical_particle_factor=False)

    reverse = ApproximateRate(rates, approx_type="nn_g", is_reverse=True,
                              use_identical_particle_factor=False)

    return forward, reverse


class ApproximateRate(Rate):
    """An approximation to a rate sequence A <--> B.  The following
    approximations are currently supported:

    * "ap_pg" : combine the A(α,γ)B and A(α,p)X(p,γ)B sequences
      into a single, effective A(α,γ)B rate.

      An example of this is combining S32(α,γ)Ar36 and
      S32(α,p)Cl35(p,γ)Ar36.

    * "nn_g" : replace the sequence A(n,γ)X(n,γ)B with an
      effective A(nn,γ)B rate.

      An example of this is combining Fe52(n,γ)Fe53(n,γ)Fe54

    * "Yp_pg" : combine A(Y,γ)B and A(Y,p)X(p,γ)B sequences into a
      single, effective A(Y,γ)B rate.  Note: the original A(Y,γ)B does
      not need to be included, in which case just A(Y,p)X(p,γ)B is
      approximated.  We also include another pathway from X connecting
      to a different nucleus C, X(p,α)C.  This affects the branching
      from X that serves as the normalization.

      Here Y is another nucleus, and we assume that mass_number(A) >=
      mass_number(Y)

      An example of this is combining O16(O16,p)P31(p,γ)S32 and
      optionally (a modified rate version of) O16(O16,γ)S32.

    * "Yp_pa" : combine the A(Y,α)B and A(Y,p)X(p,α)B sequences into a
      single, effective A(Y,α)B rate.  We also include another pathway
      from X connecting to a different nucleus C, X(p,γ)C. This affects
      the branching from X that serves as the normalization.

    This class stores all of the rates needed to implement these
    approximations.

    Parameters
    ----------
    rates : dict(str)
        A dictionary keyed by the generic form of the rate, e.g.,
        "A(a,p)X", that provides the ``Rate`` object for that rate.
        Each approximation has a different set of keys that is
        expected.
    is_reverse : bool
        Are we creating the effective A(x,y)B or B(y,x)A?
    approx_type : str
        The type of approximation to do.  Currently supported are
        "ap_pg", "nn_g", "Yp_pg", and "Yp_pa"
    use_identical_particle_factor : bool
        Usually if a rate has 2 reactants of the same type, we
        divide by 2, since the order doesn't matter.  However, for
        some approximations, like A(n,g)X(n,g)B -> A(nn,g), we
        don't want this factor, since the neutron captures are
        sequential.  This option allows the double counting
        factor to be disabled.

    """

    def __init__(self, rates, *,
                 is_reverse=False, approx_type="ap_pg",
                 use_identical_particle_factor=True):

        # this will hold all of the rates
        self.rates = rates

        # this will hold only those rates that are approximated out.  This is
        # used primarily for the RateCollection plot()
        self.hidden_rates = []

        self.is_reverse = is_reverse

        self.approx_type = approx_type

        # some approximate rates will require density or composition,
        # so when we output the python or C++ function, we will need
        # these in the argument list.
        self.rate_eval_needs_rho = False
        self.rate_eval_needs_comp = False

        if self.approx_type == "ap_pg":

            # an ap_pg approximate rate combines A(a,g)B and
            # A(a,p)X(p,g)B into a single effective rate by assuming
            # proton equilibrium.  There can be 6 or 7 rates (7 for
            # the case when we are expecting another branching from X)

            assert len(self.rates) == 6 or len(self.rates) == 7

            # make sure the keys we expect are valid in the rates dict

            try:
                # this primary rate is the forward rate that connects
                # the nucleus endpoints directly. This should be A(a,g)B
                primary_rate = self.rates["A(a,g)B"]
            except KeyError:
                print("primary rate not found")
                raise

            _assert_rate_prop(primary_rate,
                              reactants=[Nucleus("he4")], num_products=1)

            # we are going to define the product A and reactant B from
            # this reaction
            self.primary_reactant = max(primary_rate.reactants)
            self.primary_product = max(primary_rate.products)

            try:
                # the first secondary rate should be A(a,p)X, where X
                # is the intermediate nucleus
                secondary_rate_1 = self.rates["A(a,p)X"]
            except KeyError:
                print("first secondary rate not found")
                raise

            _assert_rate_prop(secondary_rate_1,
                              reactants=[self.primary_reactant, Nucleus("he4")],
                              products=[Nucleus("p")])

            # get the intermediate nucleus (X) from this rate
            self.intermediate_nucleus = max(secondary_rate_1.products)

            try:
                # now the second secondary rate show be X(p,g)B
                secondary_rate_2 = self.rates["X(p,g)B"]
            except KeyError:
                print("second secondary rate not found")
                raise

            _assert_rate_prop(secondary_rate_2,
                              reactants=[self.intermediate_nucleus, Nucleus("p")],
                              products=[self.primary_product], num_products=1)

            # now ensure that the reverse rate makes sense

            try:
                # the primary reverse rate is B(g,a)A
                primary_reverse = self.rates["B(g,a)A"]
            except KeyError:
                print("primary reverse rate not found")
                raise

            _assert_rate_prop(primary_reverse,
                              reactants=[self.primary_product],
                              products=[self.primary_reactant])

            try:
                # the first secondary reverse rate should be B(g,p)X
                secondary_reverse_1 = self.rates["B(g,p)X"]
            except KeyError:
                print("first secondary reverse rate not found")
                raise

            _assert_rate_prop(secondary_reverse_1,
                              reactants=[self.primary_product],
                              products=[self.intermediate_nucleus, Nucleus("p")])

            try:
                # the second secondary reverse rate should be X(p,a)A
                secondary_reverse_2 = self.rates["X(p,a)A"]
            except KeyError:
                print("second secondary reverse rate not found")
                raise

            _assert_rate_prop(secondary_reverse_2,
                              reactants=[self.intermediate_nucleus, Nucleus("p")],
                              products=[self.primary_reactant, Nucleus("he4")])

            # finally, check to see if we have the additional branching rate
            alternate_rate = None
            try:
                # the alternate branching rate should be X(p,Y)C
                alternate_rate = self.rates["X(p,Y)C"]
            except KeyError:
                if len(self.rates) != 6:
                    print("alternate branching rate X(p,Y)C missing")
                    raise

            if alternate_rate:
                _assert_rate_prop(alternate_rate,
                                  reactants=[self.intermediate_nucleus, Nucleus("p")],
                                  num_products=2)

            # now initialize the super class with these reactants and products

            if not self.is_reverse:
                super().__init__(reactants=[self.primary_reactant, Nucleus("he4")],
                                 products=[self.primary_product],
                                 label="approx",
                                 use_identical_particle_factor=use_identical_particle_factor)
            else:
                super().__init__(reactants=[self.primary_product],
                                 products=[self.primary_reactant, Nucleus("he4")],
                                 label="approx",
                                 use_identical_particle_factor=use_identical_particle_factor)

            self.hidden_rates = [self.rates["A(a,p)X"],
                                 self.rates["X(p,g)B"],
                                 self.rates["B(g,p)X"],
                                 self.rates["X(p,a)A"]]

        elif self.approx_type == "nn_g":

            # a nn_g approximate rate combines A(n,g)X(n,g)B into a
            # single effective rate by assuming equilibrium of X.

            assert len(self.rates) == 4

            # make sure that the pair of forward rates makes sense

            try:
                # the first forward rate should be A(n,g)X
                forward1 = self.rates["A(n,g)X"]
            except KeyError:
                print("first forward rate not found")
                raise

            _assert_rate_prop(forward1,
                              reactants=[Nucleus("n")], num_products=1)

            try:
                # the second forward rate should be X(n,g)B
                forward2 = self.rates["X(n,g)B"]
            except KeyError:
                print("second forward rate not found")
                raise

            _assert_rate_prop(forward2,
                              reactants=[Nucleus("n")], num_products=1)

            # make sure that the intermediate nucleus matches
            assert forward1.products[0] == max(forward2.reactants)

            # we are going to define the product A and reactant B from
            # these forward rates

            self.primary_reactant = max(forward1.reactants)
            self.primary_product = max(forward2.products)

            # also get the intermediate nucleus

            self.intermediate_nucleus = max(forward1.products)

            # now ensure that the reverse rates makes sense

            try:
                # the first reverse rate should be B(g,n)X
                reverse1 = self.rates["B(g,n)X"]
            except KeyError:
                print("first reverse rate not found")
                raise

            _assert_rate_prop(reverse1,
                              reactants=[self.primary_product], num_reactants=1,
                              products=[self.intermediate_nucleus, Nucleus("n")])

            try:
                # the second reverse rate should be X(g,n)A
                reverse2 = self.rates["X(g,n)A"]
            except KeyError:
                print("second reverse rate not found")
                raise

            _assert_rate_prop(reverse2,
                              reactants=[self.intermediate_nucleus], num_reactants=1,
                              products=[self.primary_reactant, Nucleus("n")])

            # now initialize the super class with these reactants and products

            # this only makes sense with the use_identical_particle_factor = False
            assert not use_identical_particle_factor

            if not self.is_reverse:
                super().__init__(reactants=[self.primary_reactant, Nucleus("n"), Nucleus("n")],
                                 products=[self.primary_product],
                                 label="approx",
                                 use_identical_particle_factor=use_identical_particle_factor)
            else:
                super().__init__(reactants=[self.primary_product],
                                 products=[self.primary_reactant, Nucleus("n"), Nucleus("n")],
                                 label="approx",
                                 use_identical_particle_factor=use_identical_particle_factor)

            self.approx = True

            # none of these rates directly appear as links in the network
            self.hidden_rates = list(self.rates.values())

            self.rate_eval_needs_rho = True
            self.rate_eval_needs_comp = True

        elif self.approx_type == "Yp_pg":

            try:
                # this is the first forward rate
                first_forward = self.rates["A(Y,p)X"]
            except KeyError:
                print("rate A(Y,p)X not found")
                raise

            _assert_rate_prop(first_forward, products=[Nucleus("p")])

            # get the primary reactant, other reactant (Y), and
            # intermediate nucleus (X).
            self.primary_reactant = max(first_forward.reactants)
            self.other_reactant = min(first_forward.reactants)
            self.intermediate_nucleus = max(first_forward.products)

            try:
                # this is the second forward rate
                second_forward = self.rates["X(p,g)B"]
            except KeyError:
                print("rate X(p,g)B not found")
                raise

            _assert_rate_prop(second_forward,
                              reactants=[self.intermediate_nucleus, Nucleus("p")],
                              num_products=1)

            # get the primary product from this rate
            self.primary_product = max(second_forward.products)

            try:
                # this is the first reverse rate
                first_reverse = self.rates["B(g,p)X"]
            except KeyError:
                print("rate B(g,p)X not found")
                raise

            _assert_rate_prop(first_reverse,
                              reactants=[self.primary_product], num_reactants=1,
                              products=[self.intermediate_nucleus, Nucleus("p")])

            try:
                # this is the second reverse rate
                second_reverse = self.rates["X(p,Y)A"]
            except KeyError:
                print("rate X(p,Y)A not found")
                raise

            _assert_rate_prop(second_reverse,
                              reactants=[self.intermediate_nucleus, Nucleus("p")],
                              products=[self.primary_reactant, self.other_reactant])

            try:
                # This is the direct rate that may optionally be present
                # it connects the nucleus endpoints
                direct_rate = self.rates["A(Y,g)B"]
            except KeyError:
                direct_rate = None

            if direct_rate:
                _assert_rate_prop(direct_rate,
                                  reactants=[self.primary_reactant, self.other_reactant],
                                  products=[self.primary_product], num_products=1)

            try:
                # This is the direct reverse rate that may be optionally present
                direct_reverse = self.rates["B(g,Y)A"]
            except KeyError:
                direct_reverse = None

            if direct_reverse:
                _assert_rate_prop(direct_reverse,
                                  reactants=[self.primary_product], num_reactants=1,
                                  products=[self.primary_reactant, self.other_reactant])

            assert (direct_rate and direct_reverse) or (not direct_rate and not direct_reverse)

            try:
                # The alternate branching rate
                alternate_branch = self.rates["X(p,a)C"]
            except KeyError:
                print("rate X(p,a)C not found")
                raise

            _assert_rate_prop(alternate_branch,
                              reactants=[self.intermediate_nucleus, Nucleus("p")],
                              products=[Nucleus("he4")])

            if not self.is_reverse:
                super().__init__(reactants=[self.primary_reactant, self.other_reactant],
                                 products=[self.primary_product],
                                 label="approx",
                                 use_identical_particle_factor=use_identical_particle_factor)
                self.hidden_rates = [self.rates["A(Y,p)X"],
                                     self.rates["X(p,g)B"],
                                     self.rates["X(p,Y)A"],
                                     self.rates["X(p,a)C"]]

            else:
                super().__init__(reactants=[self.primary_product],
                                 products=[self.primary_reactant, self.other_reactant],
                                 label="approx",
                                 use_identical_particle_factor=use_identical_particle_factor)

                self.hidden_rates = [self.rates["X(p,g)B"],
                                     self.rates["B(g,p)X"],
                                     self.rates["X(p,Y)A"],
                                     self.rates["X(p,a)C"]]

        elif self.approx_type == "Yp_pa":

            try:
                # this is the single-step forward rate
                primary_forward = self.rates["A(Y,a)B"]
            except KeyError:
                print("rate A(Y,a)B not found")
                raise

            _assert_rate_prop(primary_forward, products=[Nucleus("he4")])

            # get the primary reactant, other reactant (Y), and
            # primary product
            self.primary_reactant = max(primary_forward.reactants)
            self.other_reactant = min(primary_forward.reactants)
            self.primary_product = max(primary_forward.products)

            try:
                # this is the forward rate 1 for the 2-rate sequence
                first_forward = self.rates["A(Y,p)X"]
            except KeyError:
                print("rate A(Y,p)X not found")
                raise

            _assert_rate_prop(first_forward,
                              reactants=[self.primary_reactant, self.other_reactant],
                              products=[Nucleus("p")])

            # get the intermediate nucleus
            self.intermediate_nucleus = max(first_forward.products)

            try:
                # this is the forward rate 2 for the 2-rate sequence
                second_forward = self.rates["X(p,a)B"]
            except KeyError:
                print("rate X(p,g)B not found")
                raise

            _assert_rate_prop(second_forward,
                              reactants=[self.intermediate_nucleus, Nucleus("p")],
                              products=[Nucleus("he4"), self.primary_product])

            try:
                # this is the single-step reverse rate
                primary_reverse = rates["B(a,Y)A"]
            except KeyError:
                print("rate B(a,Y)A not found")
                raise

            _assert_rate_prop(primary_reverse,
                              reactants=[self.primary_product, Nucleus("he4")],
                              products=[self.primary_reactant, self.other_reactant])

            try:
                # this is the reverse rate 1 for the 2-rate sequence
                first_reverse = self.rates["B(a,p)X"]
            except KeyError:
                print("rate B(a,p)X not found")
                raise

            _assert_rate_prop(first_reverse,
                              reactants=[self.primary_product, Nucleus("he4")],
                              products=[self.intermediate_nucleus, Nucleus("p")])

            try:
                # this is the reverse rate 2 for the 2-rate sequence
                second_reverse = self.rates["X(p,Y)A"]
            except KeyError:
                print("rate X(p,Y)A not found")
                raise

            _assert_rate_prop(second_reverse,
                              reactants=[self.intermediate_nucleus, Nucleus("p")],
                              products=[self.primary_reactant, self.other_reactant])

            try:
                # The alternate branching rate
                alternate_branch = self.rates["X(p,g)C"]
            except KeyError:
                print("rate X(p,g)C not found")
                raise

            _assert_rate_prop(alternate_branch,
                              reactants=[self.intermediate_nucleus, Nucleus("p")],
                              num_products=1)

            if not self.is_reverse:
                super().__init__(reactants=[self.primary_reactant, self.other_reactant],
                                 products=[self.primary_product, Nucleus("he4")],
                                 label="approx",
                                 use_identical_particle_factor=use_identical_particle_factor)
                self.hidden_rates = [self.rates["A(Y,p)X"],
                                     self.rates["X(p,a)B"],
                                     self.rates["X(p,g)C"],
                                     self.rates["X(p,Y)A"]]

            else:
                super().__init__(reactants=[self.primary_product, Nucleus("he4")],
                                 products=[self.primary_reactant, self.other_reactant],
                                 label="approx",
                                 use_identical_particle_factor=use_identical_particle_factor)

                self.hidden_rates = [self.rates["B(a,p)X"],
                                     self.rates["X(p,a)B"],
                                     self.rates["X(p,Y)A"],
                                     self.rates["X(p,g)C"]]

        else:
            raise NotImplementedError(f"approximation type {self.approx_type} not supported")

        # update the Q value
        self._set_q()

    def get_child_rates(self):
        """Return a list of all of the rates that are used in this
        approximation.

        Returns
        -------
        list(Rate)

        """
        return list(self.rates.values())

    def _set_screening(self):
        # the individual rates are screened -- we don't screen the combination of them
        self.ion_screen = []
        self.screening_pairs = []

    def log_eval(self, T, *, rho=None, comp=None,
                 screen_func=None):
        """Evaluate the natural log of reaction rate for approximate rate.

        Parameters
        ----------
        T : float
            the temperature to evaluate the rate at
        rho : float
            the density to evaluate screening effects at.
        comp : float
            the composition (of type
            :py:class:`Composition <pynucastro.nucdata.composition.Composition>`)
            to evaluate screening effects with.
        screen_func : Callable
            one of the screening functions from :py:mod:`pynucastro.screening`
            -- if provided, then the rate will include screening correction.

        Returns
        -------
        float
        """

        return math.log(self.eval(T, rho=rho, comp=comp, screen_func=screen_func))

    def eval(self, T, *, rho=None, comp=None,  # pylint: disable=too-many-return-statements
             screen_func=None):
        """Evaluate the approximate rate.

        Parameters
        ----------
        T : float
            the temperature to evaluate the rate at
        rho : float
            the density to evaluate screening effects at.
        comp : float
            the composition (of type
            :py:class:`Composition <pynucastro.nucdata.composition.Composition>`)
            to evaluate screening effects with.
        screen_func : Callable
            one of the screening functions from :py:mod:`pynucastro.screening`
            -- if provided, then the rate will include screening correction.

        Returns
        -------
        float
        """

        if self.approx_type == "ap_pg":

            r_pg = self.rates["X(p,g)B"].eval(T, rho=rho, comp=comp,
                                              screen_func=screen_func)
            r_pa = self.rates["X(p,a)A"].eval(T, rho=rho, comp=comp,
                                              screen_func=screen_func)
            # there might be another branching from X that we need to
            # account for
            r_pY = 0.0
            if "X(p,Y)C" in self.rates:
                r_pY = self.rates["X(p,Y)C"].eval(T, rho=rho, comp=comp,
                                                  screen_func=screen_func)

            denom = r_pg + r_pa + r_pY

            if not self.is_reverse:  # pylint: disable=no-else-return
                # the approximate forward rate is r_ag + r_ap r_pg / (r_pg + r_pa)
                r_ag = self.rates["A(a,g)B"].eval(T, rho=rho, comp=comp,
                                                  screen_func=screen_func)
                r_ap = self.rates["A(a,p)X"].eval(T, rho=rho, comp=comp,
                                                  screen_func=screen_func)
                return r_ag + r_ap * r_pg / denom

            else:
                # the approximate reverse rate is r_ga + r_pa r_gp / (r_pg + r_pa)

                r_ga = self.rates["B(g,a)A"].eval(T, rho=rho, comp=comp,
                                                  screen_func=screen_func)
                r_gp = self.rates["B(g,p)X"].eval(T, rho=rho, comp=comp,
                                                  screen_func=screen_func)
                return r_ga + r_pa * r_gp / denom

        elif self.approx_type == "nn_g":

            # we are approximating A(n,g)X(n,g)B

            Yn = comp.get_molar()[Nucleus("n")]
            X_ng_B = self.rates["X(n,g)B"].eval(T, rho=rho, comp=comp,
                                                screen_func=screen_func)
            X_gn_A = self.rates["X(g,n)A"].eval(T, rho=rho, comp=comp,
                                                screen_func=screen_func)
            denom = rho * Yn * X_ng_B + X_gn_A

            if not self.is_reverse:  # pylint: disable=no-else-return
                # the forward rate
                A_ng_X = self.rates["A(n,g)X"].eval(T, rho=rho, comp=comp,
                                                    screen_func=screen_func)
                return A_ng_X * X_ng_B / denom

            else:
                # the reverse rate
                B_gn_X = self.rates["B(g,n)X"].eval(T, rho=rho, comp=comp,
                                                    screen_func=screen_func)
                return B_gn_X * X_gn_A / denom

        elif self.approx_type == "Yp_pg":

            # we are approximating A(Y,p)X(p,g)B with an alternate
            # branch from X, X(p,a)C, and possibly a direct path
            # between A and B, A(Y,g)B

            r_pY = self.rates["X(p,Y)A"].eval(T, rho=rho, comp=comp,
                                              screen_func=screen_func)
            r_pa = self.rates["X(p,a)C"].eval(T, rho=rho, comp=comp,
                                              screen_func=screen_func)
            r_pg = self.rates["X(p,g)B"].eval(T, rho=rho, comp=comp,
                                              screen_func=screen_func)
            denom = r_pY + r_pa + r_pg

            if not self.is_reverse:  # pylint: disable=no-else-return
                # optional direct forward rate
                r_Yg = 0.0
                try:
                    rr = self.rates["A(Y,g)B"]
                    r_Yg = rr.eval(T, rho=rho, comp=comp,
                                   screen_func=screen_func)
                except KeyError:
                    pass

                # forward rates
                r_Yp = self.rates["A(Y,p)X"].eval(T, rho=rho, comp=comp,
                                                  screen_func=screen_func)
                return r_Yg + r_Yp * r_pg / denom

            else:
                # optional direct reverse rate
                r_gY = 0.0
                try:
                    rr = self.rates["B(g,Y)A"]
                    r_gY = rr.eval(T, rho=rho, comp=comp,
                                   screen_func=screen_func)
                except KeyError:
                    pass

                # reverse
                r_gp = self.rates["B(g,p)X"].eval(T, rho=rho, comp=comp,
                                                  screen_func=screen_func)
                return r_gY + r_pY * r_gp / denom

        elif self.approx_type == "Yp_pa":

            # we are approximating A(Y,a)B + A(Y,p)X(p,a)B with an alternate
            # branch from X, X(p,g)C

            r_pY = self.rates["X(p,Y)A"].eval(T, rho=rho, comp=comp,
                                              screen_func=screen_func)
            r_pa = self.rates["X(p,a)B"].eval(T, rho=rho, comp=comp,
                                              screen_func=screen_func)
            r_pg = self.rates["X(p,g)C"].eval(T, rho=rho, comp=comp,
                                              screen_func=screen_func)
            denom = r_pY + r_pa + r_pg

            if not self.is_reverse:  # pylint: disable=no-else-return
                # direct forward rate
                r_Ya = self.rates["A(Y,a)B"].eval(T, rho=rho, comp=comp,
                                                  screen_func=screen_func)

                # forward rates
                r_Yp = self.rates["A(Y,p)X"].eval(T, rho=rho, comp=comp,
                                                  screen_func=screen_func)
                return r_Ya + r_Yp * r_pa / denom

            else:
                # direct reverse rate
                r_aY = self.rates["B(a,Y)A"].eval(T, rho=rho, comp=comp,
                                                  screen_func=screen_func)
                # reverse rates
                r_ap = self.rates["B(a,p)X"].eval(T, rho=rho, comp=comp,
                                                  screen_func=screen_func)
                return r_aY + r_pY * r_ap / denom

        raise NotImplementedError(f"approximation type {self.approx_type} not supported")

    def function_string_py(self):
        """Return a string containing the python function that
        computes the approximate rate.

        Returns
        -------
        str

        """

        if self.approx_type == "ap_pg":

            # we are approximating A(a,p)X(p,g)B
            # with the possibility of X(p,Y)C as another branching

            string = ""
            string += "@numba.njit()\n"
            string += f"def {self.fname}(rate_eval, tf):\n"

            string += f"    r_pg = rate_eval.{self.rates['X(p,g)B'].fname}\n"
            string += f"    r_pa = rate_eval.{self.rates['X(p,a)A'].fname}\n"
            if "X(p,Y)C" in self.rates:
                string += f"    r_pY = rate_eval.{self.rates['X(p,Y)C'].fname}\n"
            else:
                string += "    r_pY = 0.0\n"

            if not self.is_reverse:

                # first we need to get all of the rates that make this up
                string += f"    r_ag = rate_eval.{self.rates['A(a,g)B'].fname}\n"
                string += f"    r_ap = rate_eval.{self.rates['A(a,p)X'].fname}\n"

                # now the approximation
                string += "    rate = r_ag + r_ap * r_pg / (r_pg + r_pa + r_pY)\n"

            else:

                # first we need to get all of the rates that make this up
                string += f"    r_ga = rate_eval.{self.rates['B(g,a)A'].fname}\n"
                string += f"    r_gp = rate_eval.{self.rates['B(g,p)X'].fname}\n"

                # now the approximation
                string += "    rate = r_ga + r_pa * r_gp / (r_pg + r_pa + r_pY)\n"

            string += f"    rate_eval.{self.fname} = rate\n\n"
            return string

        if self.approx_type == "nn_g":

            # we are approximating A(n,g)X(n,g)B

            string = ""
            string += "@numba.njit()\n"
            string += f"def {self.fname}(rate_eval, tf, rho=None, Y=None):\n"

            string += "    Yn = Y[jn]\n"

            if not self.is_reverse:

                # first we need to get all of the rates that make this up
                string += f"    r1_ng = rate_eval.{self.rates['A(n,g)X'].fname}\n"
                string += f"    r2_ng = rate_eval.{self.rates['X(n,g)B'].fname}\n"
                string += f"    r1_gn = rate_eval.{self.rates['X(g,n)A'].fname}\n"

                # now the approximation
                string += "    rate = r1_ng * r2_ng / (rho * Yn * r2_ng + r1_gn)\n"

            else:

                # first we need to get all of the rates that make this up
                string += f"    r1_gn = rate_eval.{self.rates['X(g,n)A'].fname}\n"
                string += f"    r2_gn = rate_eval.{self.rates['B(g,n)X'].fname}\n"
                string += f"    r2_ng = rate_eval.{self.rates['X(n,g)B'].fname}\n"

                # now the approximation
                string += "    rate = r1_gn * r2_gn / (rho * Yn * r2_ng + r1_gn)\n"

            string += f"    rate_eval.{self.fname} = rate\n\n"
            return string

        if self.approx_type == "Yp_pg":

            # we are approximating A(Y,p)X(p,g)B with an alternate
            # branch from X, X(p,a)C, and possibly a direct path
            # between A and B, A(Y,g)B

            string = ""
            string += "@numba.njit()\n"
            string += f"def {self.fname}(rate_eval, tf):\n"

            string += f"    r_pY = rate_eval.{self.rates['X(p,Y)A'].fname}\n"
            string += f"    r_pa = rate_eval.{self.rates['X(p,a)C'].fname}\n"
            string += f"    r_pg = rate_eval.{self.rates['X(p,g)B'].fname}\n"
            string += "    denom = r_pY + r_pa + r_pg\n"

            if not self.is_reverse:

                # first we need to get all of the rates that make this up
                if "A(Y,g)B" in self.rates:
                    string += f"    r_Yg = rate_eval.{self.rates['A(Y,g)B'].fname}\n"
                else:
                    string += "    r_Yg = 0.0\n"

                string += f"    r_Yp = rate_eval.{self.rates['A(Y,p)X'].fname}\n"

                # now the approximation
                string += "    rate = r_Yg + r_Yp * r_pg / denom\n"

            else:

                # first we need to get all of the rates that make this up
                if "B(g,Y)A" in self.rates:
                    string += f"    r_gY = rate_eval.{self.rates['B(g,Y)A'].fname}\n"
                else:
                    string += "    r_gY = 0.0\n"

                string += f"    r_gp = rate_eval.{self.rates['B(g,p)X'].fname}\n"

                # now the approximation
                string += "    rate = r_gY + r_pY * r_gp / denom\n"

            string += f"    rate_eval.{self.fname} = rate\n\n"
            return string

        if self.approx_type == "Yp_pa":

            # we are approximating A(Y,a)B + A(Y,p)X(p,a)B with an alternate
            # branch from X, X(p,g)C

            string = ""
            string += "@numba.njit()\n"
            string += f"def {self.fname}(rate_eval, tf):\n"

            string += f"    r_pY = rate_eval.{self.rates['X(p,Y)A'].fname}\n"
            string += f"    r_pa = rate_eval.{self.rates['X(p,a)B'].fname}\n"
            string += f"    r_pg = rate_eval.{self.rates['X(p,g)C'].fname}\n"
            string += "    denom = r_pY + r_pa + r_pg\n"

            if not self.is_reverse:

                # first we need to get all of the rates that make this up
                string += f"    r_Ya = rate_eval.{self.rates['A(Y,a)B'].fname}\n"
                string += f"    r_Yp = rate_eval.{self.rates['A(Y,p)X'].fname}\n"

                # now the approximation
                string += "    rate = r_Ya + r_Yp * r_pa / denom\n"

            else:

                # first we need to get all of the rates that make this up
                string += f"    r_aY = rate_eval.{self.rates['B(a,Y)A'].fname}\n"
                string += f"    r_ap = rate_eval.{self.rates['B(a,p)X'].fname}\n"

                # now the approximation
                string += "    rate = r_aY + r_pY * r_ap / denom\n"

            string += f"    rate_eval.{self.fname} = rate\n\n"
            return string

        raise NotImplementedError("don't know how to work with this approximation")

    def function_string_cxx(self, dtype="double", specifiers="inline",
                            leave_open=False, extra_args=None):
        """Return a string containing the C++ function that computes
        the approximate rate

        Parameters
        ----------
        dtype : str
            The C++ datatype to use for all declarations
        specifiers : str
            C++ specifiers to add before each function declaration
            (i.e. "inline")
        leave_open : bool
            If ``true``, then we leave the function unclosed (no "}"
            at the end).  This can allow additional functions to add
            to this output.
        extra_args : list(str)
            A list of strings representing additional arguments that
            should be appended to the argument list when defining the
            function interface.

        Returns
        -------
        str

        """

        if extra_args is None:
            extra_args = ()

        if dtype == "amrex::Real":
            array_type = "amrex::Array1D"
        else:
            array_type = "Array1D"

        if self.approx_type == "ap_pg":

            args = ["const T& rate_eval", f"{dtype}& rate", f"{dtype}& drate_dT", *extra_args]
            fstring = ""
            fstring = "template <typename T>\n"
            fstring += f"{specifiers}\n"
            fstring += f"void rate_{self.fname}({', '.join(args)}) {{\n\n"

            fstring += f"    {dtype} r_pg = rate_eval.screened_rates(k_{self.rates['X(p,g)B'].fname});\n"
            fstring += f"    {dtype} r_pa = rate_eval.screened_rates(k_{self.rates['X(p,a)A'].fname});\n"
            if "X(p,Y)C" in self.rates:
                fstring += f"    {dtype} r_pY = rate_eval.screened_rates(k_{self.rates['X(p,Y)C'].fname});\n"
            else:
                fstring += f"    {dtype} r_pY = 0.0_rt;\n"

            fstring += f"    {dtype} dd = 1.0_rt / (r_pg + r_pa + r_pY);\n"

            if not self.is_reverse:

                # first we need to get all of the rates that make this up
                fstring += f"    {dtype} r_ag = rate_eval.screened_rates(k_{self.rates['A(a,g)B'].fname});\n"
                fstring += f"    {dtype} r_ap = rate_eval.screened_rates(k_{self.rates['A(a,p)X'].fname});\n"

                # now the approximation
                fstring += "    rate = r_ag + r_ap * r_pg * dd;\n"
                fstring += "    if constexpr (std::is_same_v<T, rate_derivs_t>) {\n"
                fstring += f"        {dtype} drdT_ag = rate_eval.dscreened_rates_dT(k_{self.rates['A(a,g)B'].fname});\n"
                fstring += f"        {dtype} drdT_ap = rate_eval.dscreened_rates_dT(k_{self.rates['A(a,p)X'].fname});\n"
                fstring += f"        {dtype} drdT_pg = rate_eval.dscreened_rates_dT(k_{self.rates['X(p,g)B'].fname});\n"
                fstring += f"        {dtype} drdT_pa = rate_eval.dscreened_rates_dT(k_{self.rates['X(p,a)A'].fname});\n"
                if "X(p,Y)C" in self.rates:
                    fstring += f"        {dtype} drdT_pY = rate_eval.dscreened_rates_dT(k_{self.rates['X(p,Y)C'].fname});\n"
                else:
                    fstring += f"        {dtype} drdT_pY = 0.0_rt;\n"

                fstring += "        drate_dT = drdT_ag + drdT_ap * r_pg * dd + r_ap * drdT_pg * dd - r_ap * r_pg * dd * dd * (drdT_pg + drdT_pa + drdT_pY);\n"
                fstring += "    }\n"
            else:

                # first we need to get all of the rates that make this up
                fstring += f"    {dtype} r_ga = rate_eval.screened_rates(k_{self.rates['B(g,a)A'].fname});\n"
                fstring += f"    {dtype} r_gp = rate_eval.screened_rates(k_{self.rates['B(g,p)X'].fname});\n"

                # now the approximation
                fstring += "    rate = r_ga + r_gp * r_pa * dd;\n"
                fstring += "    if constexpr (std::is_same_v<T, rate_derivs_t>) {\n"
                fstring += f"        {dtype} drdT_ga = rate_eval.dscreened_rates_dT(k_{self.rates['B(g,a)A'].fname});\n"
                fstring += f"        {dtype} drdT_pa = rate_eval.dscreened_rates_dT(k_{self.rates['X(p,a)A'].fname});\n"
                fstring += f"        {dtype} drdT_gp = rate_eval.dscreened_rates_dT(k_{self.rates['B(g,p)X'].fname});\n"
                fstring += f"        {dtype} drdT_pg = rate_eval.dscreened_rates_dT(k_{self.rates['X(p,g)B'].fname});\n"
                if "X(p,Y)C" in self.rates:
                    fstring += f"        {dtype} drdT_pY = rate_eval.dscreened_rates_dT(k_{self.rates['X(p,Y)C'].fname});\n"
                else:
                    fstring += f"        {dtype} drdT_pY = 0.0_rt;\n"

                fstring += "        drate_dT = drdT_ga + drdT_gp * r_pa * dd + r_gp * drdT_pa * dd - r_gp * r_pa * dd * dd * (drdT_pg + drdT_pa + drdT_pY);\n"
                fstring += "    }\n"

            if not leave_open:
                fstring += "}\n\n"

            return fstring

        if self.approx_type == "nn_g":

            args = ["const T& rate_eval", f"const {dtype} rho", f"const {array_type}<{dtype}, 1, NumSpec>& Y",
                    f"{dtype}& rate", f"{dtype}& drate_dT", *extra_args]
            fstring = ""
            fstring = "template <typename T>\n"
            fstring += f"{specifiers}\n"
            fstring += f"void rate_{self.fname}({', '.join(args)}) {{\n\n"
            fstring += f"    {dtype} Yn = Y(N);\n"

            if not self.is_reverse:

                # first we need to get all of the rates that make this up
                fstring += f"    {dtype} r1_ng = rate_eval.screened_rates(k_{self.rates['A(n,g)X'].fname});\n"
                fstring += f"    {dtype} r2_ng = rate_eval.screened_rates(k_{self.rates['X(n,g)B'].fname});\n"
                fstring += f"    {dtype} r1_gn = rate_eval.screened_rates(k_{self.rates['X(g,n)A'].fname});\n"

                # now the approximation
                fstring += f"    {dtype} dd = 1.0_rt / (rho * Yn * r2_ng + r1_gn);\n"
                fstring += "    rate = r1_ng * r2_ng * dd;\n"
                fstring += "    if constexpr (std::is_same_v<T, rate_derivs_t>) {\n"
                fstring += f"        {dtype} dr1dT_ng = rate_eval.dscreened_rates_dT(k_{self.rates['A(n,g)X'].fname});\n"
                fstring += f"        {dtype} dr2dT_ng = rate_eval.dscreened_rates_dT(k_{self.rates['X(n,g)B'].fname});\n"
                fstring += f"        {dtype} dr1dT_gn = rate_eval.dscreened_rates_dT(k_{self.rates['X(g,n)A'].fname});\n"
                fstring += "        drate_dT = dr1dT_ng * r2_ng * dd + r1_ng * dr2dT_ng * dd - r1_ng * r2_ng * dd * dd * (rho * Yn * dr2dT_ng + dr1dT_gn);\n"
                fstring += "    }\n"
            else:

                # first we need to get all of the rates that make this up
                fstring += f"    {dtype} r1_gn = rate_eval.screened_rates(k_{self.rates['X(g,n)A'].fname});\n"
                fstring += f"    {dtype} r2_gn = rate_eval.screened_rates(k_{self.rates['B(g,n)X'].fname});\n"
                fstring += f"    {dtype} r2_ng = rate_eval.screened_rates(k_{self.rates['X(n,g)B'].fname});\n"

                # now the approximation
                fstring += f"    {dtype} dd = 1.0_rt / (rho * Yn * r2_ng + r1_gn);\n"
                fstring += "    rate = r1_gn * r2_gn * dd;\n"
                fstring += "    if constexpr (std::is_same_v<T, rate_derivs_t>) {\n"
                fstring += f"        {dtype} dr1dT_gn = rate_eval.dscreened_rates_dT(k_{self.rates['X(g,n)A'].fname});\n"
                fstring += f"        {dtype} dr2dT_gn = rate_eval.dscreened_rates_dT(k_{self.rates['B(g,n)X'].fname});\n"
                fstring += f"        {dtype} dr2dT_ng = rate_eval.dscreened_rates_dT(k_{self.rates['X(n,g)B'].fname});\n"
                fstring += "        drate_dT = dr1dT_gn * r2_gn * dd + r1_gn * dr2dT_gn * dd - r1_gn * r2_gn * dd * dd * (rho * Yn * dr2dT_ng + dr1dT_gn);\n"
                fstring += "    }\n"

            if not leave_open:
                fstring += "}\n\n"

            return fstring

        if self.approx_type == "Yp_pg":

            # we are approximating A(Y,p)X(p,g)B with an alternate
            # branch from X, X(p,a)C, and possibly a direct path
            # between A and B, A(Y,g)B

            args = ["const T& rate_eval", f"{dtype}& rate", f"{dtype}& drate_dT", *extra_args]
            fstring = ""
            fstring = "template <typename T>\n"
            fstring += f"{specifiers}\n"
            fstring += f"void rate_{self.fname}({', '.join(args)}) {{\n\n"

            fstring += f"    {dtype} r_pY = rate_eval.screened_rates(k_{self.rates['X(p,Y)A'].fname});\n"
            fstring += f"    {dtype} r_pa = rate_eval.screened_rates(k_{self.rates['X(p,a)C'].fname});\n"
            fstring += f"    {dtype} r_pg = rate_eval.screened_rates(k_{self.rates['X(p,g)B'].fname});\n"

            fstring += f"    {dtype} dd = 1.0_rt / (r_pY + r_pa + r_pg);\n"

            if not self.is_reverse:

                # first we need to get all of the rates that make this up
                if "A(Y,g)B" in self.rates:
                    fstring += f"    {dtype} r_Yg = rate_eval.screened_rates(k_{self.rates['A(Y,g)B'].fname});\n"
                else:
                    fstring += f"    {dtype} r_Yg = 0.0_rt;\n"

                fstring += f"    {dtype} r_Yp = rate_eval.screened_rates(k_{self.rates['A(Y,p)X'].fname});\n"

                # now the approximation
                fstring += "    rate = r_Yg + r_Yp * r_pg * dd;\n"
                fstring += "    if constexpr (std::is_same_v<T, rate_derivs_t>) {\n"
                fstring += f"        {dtype} drdT_pY = rate_eval.dscreened_rates_dT(k_{self.rates['X(p,Y)A'].fname});\n"
                fstring += f"        {dtype} drdT_pa = rate_eval.dscreened_rates_dT(k_{self.rates['X(p,a)C'].fname});\n"
                fstring += f"        {dtype} drdT_pg = rate_eval.dscreened_rates_dT(k_{self.rates['X(p,g)B'].fname});\n"

                if "A(Y,g)B" in self.rates:
                    fstring += f"        {dtype} drdT_Yg = rate_eval.dscreened_rates_dT(k_{self.rates['A(Y,g)B'].fname});\n"
                else:
                    fstring += f"        {dtype} drdT_Yg = 0.0_rt;\n"
                fstring += f"        {dtype} drdT_Yp = rate_eval.dscreened_rates_dT(k_{self.rates['A(Y,p)X'].fname});\n"

                fstring += "        drate_dT = drdT_Yg + drdT_Yp * r_pg * dd + r_Yp * drdT_pg * dd - r_Yp * r_pg * dd * dd * (drdT_pY + drdT_pa + drdT_pg);\n"
                fstring += "    }\n"
            else:

                # first we need to get all of the rates that make this up
                if "B(g,Y)A" in self.rates:
                    fstring += f"    {dtype} r_gY = rate_eval.screened_rates(k_{self.rates['B(g,Y)A'].fname});\n"
                else:
                    fstring += f"    {dtype} r_gY = 0.0_rt;\n"

                fstring += f"    {dtype} r_gp = rate_eval.screened_rates(k_{self.rates['B(g,p)X'].fname});\n"

                # now the approximation
                fstring += "    rate = r_gY + r_pY * r_gp * dd;\n"
                fstring += "    if constexpr (std::is_same_v<T, rate_derivs_t>) {\n"
                fstring += f"        {dtype} drdT_pY = rate_eval.dscreened_rates_dT(k_{self.rates['X(p,Y)A'].fname});\n"
                fstring += f"        {dtype} drdT_pa = rate_eval.dscreened_rates_dT(k_{self.rates['X(p,a)C'].fname});\n"
                fstring += f"        {dtype} drdT_pg = rate_eval.dscreened_rates_dT(k_{self.rates['X(p,g)B'].fname});\n"
                if "B(g,Y)A" in self.rates:
                    fstring += f"        {dtype} drdT_gY = rate_eval.screened_rates_dT(k_{self.rates['B(g,Y)A'].fname});\n"
                else:
                    fstring += f"        {dtype} drdT_gY = 0.0_rt;\n"
                fstring += f"        {dtype} drdT_gp = rate_eval.dscreened_rates_dT(k_{self.rates['B(g,p)X'].fname});\n"

                fstring += "        drate_dT = drdT_gY + drdT_pY * r_gp * dd + r_pY * drdT_gp * dd - r_pY * r_gp * dd * dd * (drdT_pY + drdT_pa + drdT_pg);\n"
                fstring += "    }\n"

            if not leave_open:
                fstring += "}\n\n"

            return fstring

        if self.approx_type == "Yp_pa":
            raise NotImplementedError("we haven't implemented Yp_pa in C++ yet")

        raise NotImplementedError("don't know how to work with this approximation")
