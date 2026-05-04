"""Classes and methods for describing approximations to rates.  Often
this involves using a group of existing rates to enforce an
equilibrium through a nucleus.

"""

import math

from pynucastro.nucdata import Nucleus
from pynucastro.rates.rate import Rate


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
    """An approximation to a rate sequence.  Two approximations are
    currently supported:

    * "ap_pg" : combine the A(a, g)B and A(a, p)X(p, g)B sequences
      into a single, effective A(a, g)B rate.

    * "nn_g" : replace the sequence A(n, g)X(n, g)B with an
      effective A(nn, g)B rate.

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
        Are we creating the effective A(x,y)B or B(y, x)A?
    approx_type : str
        The type of approximation to do.  Currently supported are
        "ap_pg" and "nn_g"
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
            # proton equilibrium.

            assert len(self.rates) == 6

            # make sure the keys we expect are valid in the rates dict

            try:
                # this primary rate is the forward rate that connects
                # the nucleus endpoints directly
                primary_rate = self.rates["A(a,g)B"]

                # make sure that the primary forward rate makes sense
                # this should be A(a,g)B
                assert Nucleus("he4") in primary_rate.reactants
                assert len(primary_rate.products) == 1
            except KeyError:
                print("primary rate not found")
                raise

            # we are going to define the product A and reactant B from
            # this reaction
            self.primary_reactant = max(primary_rate.reactants)
            self.primary_product = max(primary_rate.products)

            try:
                # the first secondary rate should be A(a,p)X, where X
                # is the intermediate nucleus
                secondary_rate_1 = self.rates["A(a,p)X"]

                assert self.primary_reactant in secondary_rate_1.reactants
                assert Nucleus("he4") in secondary_rate_1.reactants
                assert Nucleus("p") in secondary_rate_1.products
            except KeyError:
                print("first secondary rate not found")
                raise

            # get the intermediate nucleus (X) from this rate
            self.intermediate_nucleus = max(secondary_rate_1.products)

            try:
                # now the second secondary rate show be X(p,g)B
                secondary_rate_2 = self.rates["X(p,g)B"]

                assert self.intermediate_nucleus in secondary_rate_2.reactants
                assert Nucleus("p") in secondary_rate_2.reactants
                assert self.primary_product in secondary_rate_2.products
            except KeyError:
                print("second secondary rate not found")
                raise

            # now ensure that the reverse rate makes sense

            try:
                # the primary reverse rate is B(g,a)A
                primary_reverse = self.rates["B(g,a)A"]

                assert self.primary_product in primary_reverse.reactants
                assert self.primary_reactant in primary_reverse.products
            except KeyError:
                print("primary reverse rate not found")
                raise

            try:
                # the first secondary reverse rate should be B(g,p)X
                secondary_reverse_1 = self.rates["B(g,p)X"]

                assert self.primary_product in secondary_reverse_1.reactants
                assert self.intermediate_nucleus in secondary_reverse_1.products
                assert Nucleus("p") in secondary_reverse_1.products
            except KeyError:
                print("first secondary reverse rate not found")
                raise

            try:
                # the second secondary reverse rate should be X(p,a)A
                secondary_reverse_2 = self.rates["X(p,a)A"]

                assert self.intermediate_nucleus in secondary_reverse_2.reactants
                assert Nucleus("p") in secondary_reverse_2.reactants
                assert self.primary_reactant in secondary_reverse_2.products
                assert Nucleus("he4") in secondary_reverse_2.products
            except KeyError:
                print("second secondary reverse rate not found")
                raise

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
                assert Nucleus("n") in forward1.reactants
                assert len(forward1.products) == 1
            except KeyError:
                print("first forward rate not found")
                raise

            try:
                # the second forward rate should be X(n,g)B
                forward2 = self.rates["X(n,g)B"]
                assert Nucleus("n") in forward2.reactants
                assert len(forward2.products) == 1
            except KeyError:
                print("second forward rate not found")
                raise

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
                assert self.primary_product in reverse1.reactants
                assert len(reverse1.reactants) == 1
                assert self.intermediate_nucleus in reverse1.products
                assert Nucleus("n") in reverse1.products
            except KeyError:
                print("first reverse rate not found")
                raise

            try:
                # the second reverse rate should be X(g,n)A
                reverse2 = self.rates["X(g,n)A"]
                assert self.intermediate_nucleus in reverse2.reactants
                assert len(reverse2.reactants) == 1
                assert self.primary_reactant in reverse2.products
                assert Nucleus("n") in reverse2.products
            except KeyError:
                print("second reverse rate not found")
                raise

            # now initialize the super class with these reactants and products

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

    def eval(self, T, *, rho=None, comp=None,
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
            if not self.is_reverse:  # pylint: disable=no-else-return
                # the approximate forward rate is r_ag + r_ap r_pg / (r_pg + r_pa)
                r_ag = self.rates["A(a,g)B"].eval(T, rho=rho, comp=comp,
                                                  screen_func=screen_func)
                r_ap = self.rates["A(a,p)X"].eval(T, rho=rho, comp=comp,
                                                  screen_func=screen_func)
                r_pg = self.rates["X(p,g)B"].eval(T, rho=rho, comp=comp,
                                                  screen_func=screen_func)

                r_pa = self.rates["X(p,a)A"].eval(T, rho=rho, comp=comp,
                                                  screen_func=screen_func)

                return r_ag + r_ap * r_pg / (r_pg + r_pa)

            else:
                # the approximate reverse rate is r_ga + r_pa r_gp / (r_pg + r_pa)

                r_ga = self.rates["B(g,a)A"].eval(T, rho=rho, comp=comp,
                                                  screen_func=screen_func)
                r_gp = self.rates["B(g,p)X"].eval(T, rho=rho, comp=comp,
                                                  screen_func=screen_func)
                r_pa = self.rates["X(p,a)A"].eval(T, rho=rho, comp=comp,
                                                  screen_func=screen_func)

                r_pg = self.rates["X(p,g)B"].eval(T, rho=rho, comp=comp,
                                                  screen_func=screen_func)

                return r_ga + r_pa * r_gp / (r_pg + r_pa)

        elif self.approx_type == "nn_g":

            # we are approximating A(n,g)X(n,g)B

            Yn = comp.get_molar()[Nucleus("n")]

            if not self.is_reverse:  # pylint: disable=no-else-return
                # the forward rate
                A_ng_X = self.rates["A(n,g)X"].eval(T, rho=rho, comp=comp,
                                                    screen_func=screen_func)
                X_ng_B = self.rates["X(n,g)B"].eval(T, rho=rho, comp=comp,
                                                    screen_func=screen_func)

                X_gn_A = self.rates["X(g,n)A"].eval(T, rho=rho, comp=comp,
                                                    screen_func=screen_func)

                return A_ng_X * X_ng_B / (rho * Yn * X_ng_B + X_gn_A)

            else:
                # the reverse rate
                B_gn_X = self.rates["B(g,n)X"].eval(T, rho=rho, comp=comp,
                                                    screen_func=screen_func)
                X_gn_A = self.rates["X(g,n)A"].eval(T, rho=rho, comp=comp,
                                                    screen_func=screen_func)

                X_ng_B = self.rates["X(n,g)B"].eval(T, rho=rho, comp=comp,
                                                    screen_func=screen_func)

                return B_gn_X * X_gn_A / (rho * Yn * X_ng_B + X_gn_A)

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

            string = ""
            string += "@numba.njit()\n"
            string += f"def {self.fname}(rate_eval, tf):\n"

            if not self.is_reverse:

                # first we need to get all of the rates that make this up
                string += f"    r_ag = rate_eval.{self.rates['A(a,g)B'].fname}\n"
                string += f"    r_ap = rate_eval.{self.rates['A(a,p)X'].fname}\n"
                string += f"    r_pg = rate_eval.{self.rates['X(p,g)B'].fname}\n"
                string += f"    r_pa = rate_eval.{self.rates['X(p,a)A'].fname}\n"

                # now the approximation
                string += "    rate = r_ag + r_ap * r_pg / (r_pg + r_pa)\n"

            else:

                # first we need to get all of the rates that make this up
                string += f"    r_ga = rate_eval.{self.rates['B(g,a)A'].fname}\n"
                string += f"    r_pa = rate_eval.{self.rates['X(p,a)A'].fname}\n"
                string += f"    r_gp = rate_eval.{self.rates['B(g,p)X'].fname}\n"
                string += f"    r_pg = rate_eval.{self.rates['X(p,g)B'].fname}\n"

                # now the approximation
                string += "    rate = r_ga + r_pa * r_gp / (r_pg + r_pa)\n"

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

            if not self.is_reverse:

                # first we need to get all of the rates that make this up
                fstring += f"    {dtype} r_ag = rate_eval.screened_rates(k_{self.rates['A(a,g)B'].fname});\n"
                fstring += f"    {dtype} r_ap = rate_eval.screened_rates(k_{self.rates['A(a,p)X'].fname});\n"
                fstring += f"    {dtype} r_pg = rate_eval.screened_rates(k_{self.rates['X(p,g)B'].fname});\n"
                fstring += f"    {dtype} r_pa = rate_eval.screened_rates(k_{self.rates['X(p,a)A'].fname});\n"

                # now the approximation
                fstring += f"    {dtype} dd = 1.0_rt / (r_pg + r_pa);\n"
                fstring += "    rate = r_ag + r_ap * r_pg * dd;\n"
                fstring += "    if constexpr (std::is_same_v<T, rate_derivs_t>) {\n"
                fstring += f"        {dtype} drdT_ag = rate_eval.dscreened_rates_dT(k_{self.rates['A(a,g)B'].fname});\n"
                fstring += f"        {dtype} drdT_ap = rate_eval.dscreened_rates_dT(k_{self.rates['A(a,p)X'].fname});\n"
                fstring += f"        {dtype} drdT_pg = rate_eval.dscreened_rates_dT(k_{self.rates['X(p,g)B'].fname});\n"
                fstring += f"        {dtype} drdT_pa = rate_eval.dscreened_rates_dT(k_{self.rates['X(p,a)A'].fname});\n"
                fstring += "        drate_dT = drdT_ag + drdT_ap * r_pg * dd + r_ap * drdT_pg * dd - r_ap * r_pg * dd * dd * (drdT_pg + drdT_pa);\n"
                fstring += "    }\n"
            else:

                # first we need to get all of the rates that make this up
                fstring += f"    {dtype} r_ga = rate_eval.screened_rates(k_{self.rates['B(g,a)A'].fname});\n"
                fstring += f"    {dtype} r_pa = rate_eval.screened_rates(k_{self.rates['X(p,a)A'].fname});\n"
                fstring += f"    {dtype} r_gp = rate_eval.screened_rates(k_{self.rates['B(g,p)X'].fname});\n"
                fstring += f"    {dtype} r_pg = rate_eval.screened_rates(k_{self.rates['X(p,g)B'].fname});\n"

                # now the approximation
                fstring += f"    {dtype} dd = 1.0_rt / (r_pg + r_pa);\n"
                fstring += "    rate = r_ga + r_gp * r_pa * dd;\n"
                fstring += "    if constexpr (std::is_same_v<T, rate_derivs_t>) {\n"
                fstring += f"        {dtype} drdT_ga = rate_eval.dscreened_rates_dT(k_{self.rates['B(g,a)A'].fname});\n"
                fstring += f"        {dtype} drdT_pa = rate_eval.dscreened_rates_dT(k_{self.rates['X(p,a)A'].fname});\n"
                fstring += f"        {dtype} drdT_gp = rate_eval.dscreened_rates_dT(k_{self.rates['B(g,p)X'].fname});\n"
                fstring += f"        {dtype} drdT_pg = rate_eval.dscreened_rates_dT(k_{self.rates['X(p,g)B'].fname});\n"
                fstring += "        drate_dT = drdT_ga + drdT_gp * r_pa * dd + r_gp * drdT_pa * dd - r_gp * r_pa * dd * dd * (drdT_pg + drdT_pa);\n"
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

        raise NotImplementedError("don't know how to work with this approximation")
