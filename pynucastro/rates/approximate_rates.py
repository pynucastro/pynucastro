from pynucastro.nucdata import Nucleus
from pynucastro.rates.rate import Rate


def create_double_neutron_capture(lib, reactant, product):
    """a helper function that will return an ApproximateRate object
    for the "nn_g" approximation"""

    intermediate = Nucleus(f"{reactant.el}{reactant.A+1}")

    forward_1 = lib.get_rate_by_name(f"{reactant.raw}(n,){intermediate.raw}")
    forward_2 = lib.get_rate_by_name(f"{intermediate.raw}(n,){product.raw}")

    reverse_1 = lib.get_rate_by_name(f"{product.raw}(,n){intermediate.raw}")
    reverse_2 = lib.get_rate_by_name(f"{intermediate.raw}(,n){reactant.raw}")

    forward = ApproximateRate(None, [forward_1, forward_2],
                              None, [reverse_1, reverse_2],
                              approx_type="nn_g",
                              use_identical_particle_factor=False)

    reverse = ApproximateRate(None, [forward_1, forward_2],
                              None, [reverse_1, reverse_2],
                              approx_type="nn_g", is_reverse=True,
                              use_identical_particle_factor=False)

    return forward, reverse


class ApproximateRate(Rate):

    def __init__(self, primary_rate, secondary_rates,
                 primary_reverse, secondary_reverse, *,
                 is_reverse=False, approx_type="ap_pg",
                 use_identical_particle_factor=True):
        """the primary rate has the same reactants and products as the
        final approximate rate would have.  It can be None.  The
        secondary rates are ordered such that together they would give
        the same sequence

        """

        self.primary_rate = primary_rate
        self.secondary_rates = secondary_rates

        self.primary_reverse = primary_reverse
        self.secondary_reverse = secondary_reverse

        self.is_reverse = is_reverse

        self.approx_type = approx_type

        if self.approx_type == "ap_pg":

            # an ap_pg approximate rate combines A(a,g)B and A(a,p)X(p,g)B into a
            # single effective rate by assuming proton equilibrium.

            assert len(self.secondary_rates) == 2

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
                                 label="approx",
                                 use_identical_particle_factor=use_identical_particle_factor)
            else:
                super().__init__(reactants=[self.primary_product],
                                 products=[self.primary_reactant, Nucleus("he4")],
                                 label="approx",
                                 use_identical_particle_factor=use_identical_particle_factor)

            self.chapter = "a"

        elif self.approx_type == "nn_g":

            # a nn_g approximate rate combines A(n,g)X(n,g)B into a
            # single effective rate by assuming neutron equilibrium.

            assert self.primary_rate is None
            assert len(self.secondary_rates) == 2

            # make sure that the pair of forward secondary makes sense

            for rr in self.secondary_rates:
                assert Nucleus("n") in rr.reactants and len(rr.products) == 1

            # make sure that the intermediate nucleus matches
            assert self.secondary_rates[0].products[0] == max(self.secondary_rates[1].reactants)

            # we are going to define the product A and reactant B from
            # these forward secondary rates

            self.primary_reactant = max(self.secondary_rates[0].reactants)
            self.primary_product = max(self.secondary_rates[1].products)

            # the intermediate nucleus is not in our network, so make it
            # dummy

            self.intermediate_nucleus = max(self.secondary_rates[0].products)
            #self.intermediate_nucleus.dummy = True

            # now ensure that the reverse rates makes sense

            assert self.primary_reverse is None
            assert len(self.secondary_reverse) == 2

            for rr in self.secondary_reverse:
                assert len(rr.reactants) == 1

            # now the first secondary reverse rate should be B(g,n)X

            assert (self.primary_product in self.secondary_reverse[0].reactants and
                    self.intermediate_nucleus in secondary_reverse[0].products and
                    Nucleus("n") in secondary_reverse[0].products)

            # and the second secondary reverse rate should be X(g,n)A

            assert (self.intermediate_nucleus in self.secondary_reverse[1].reactants and
                    self.primary_reactant in self.secondary_reverse[1].products and
                    Nucleus("n") in self.secondary_reverse[1].products)

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

            self.chapter = "a"

        else:
            raise NotImplementedError(f"approximation type {self.approx_type} not supported")

        # update the Q value
        self._set_q()

    def get_child_rates(self):
        """return a list of all of the rates that are used in this approximation"""
        tlist = []
        if self.primary_rate:
            tlist += [self.primary_rate]
        tlist += self.secondary_rates
        if self.primary_reverse:
            tlist += [self.primary_reverse]
        tlist += self.secondary_reverse
        return tlist

    def _set_screening(self):
        # the individual rates are screened -- we don't screen the combination of them
        pass

    def eval(self, T, *, rho=None, comp=None):
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

        elif self.approx_type == "nn_g":

            # we are approximating A(n,g)X(n,g)B

            Yn = comp.get_molar()[Nucleus("n")]

            if not self.is_reverse:  # pylint: disable=no-else-return
                # the forward rate
                A_ng_X = self.secondary_rates[0].eval(T)  # A(n,g)X
                X_ng_B = self.secondary_rates[1].eval(T)  # X(n,g)B

                X_gn_A = self.secondary_reverse[1].eval(T)  # X(g,n)A

                return A_ng_X * X_ng_B / (rho * Yn * X_ng_B + X_gn_A)

            else:
                # the reverse rate
                B_gn_X = self.secondary_reverse[0].eval(T)
                X_gn_A = self.secondary_reverse[1].eval(T)

                X_ng_B = self.secondary_rates[1].eval(T)

                return B_gn_X * X_gn_A / (rho * Yn * X_ng_B + X_gn_A)

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

    def function_string_cxx(self, dtype="double", specifiers="inline", leave_open=False, extra_args=()):
        """
        Return a string containing C++ function that computes the
        approximate rate
        """

        if self.approx_type != "ap_pg":
            raise NotImplementedError("don't know how to work with this approximation")

        args = ["const T& rate_eval", f"{dtype}& rate", f"{dtype}& drate_dT", *extra_args]
        fstring = ""
        fstring = "template <typename T>\n"
        fstring += f"{specifiers}\n"
        fstring += f"void rate_{self.cname()}({', '.join(args)}) {{\n\n"

        if not self.is_reverse:

            # first we need to get all of the rates that make this up
            fstring += f"    {dtype} r_ag = rate_eval.screened_rates(k_{self.primary_rate.cname()});\n"
            fstring += f"    {dtype} r_ap = rate_eval.screened_rates(k_{self.secondary_rates[0].cname()});\n"
            fstring += f"    {dtype} r_pg = rate_eval.screened_rates(k_{self.secondary_rates[1].cname()});\n"
            fstring += f"    {dtype} r_pa = rate_eval.screened_rates(k_{self.secondary_reverse[1].cname()});\n"

            # now the approximation
            fstring += f"    {dtype} dd = 1.0_rt / (r_pg + r_pa);\n"
            fstring += "    rate = r_ag + r_ap * r_pg * dd;\n"
            fstring += "    if constexpr (std::is_same_v<T, rate_derivs_t>) {\n"
            fstring += f"        {dtype} drdT_ag = rate_eval.dscreened_rates_dT(k_{self.primary_rate.cname()});\n"
            fstring += f"        {dtype} drdT_ap = rate_eval.dscreened_rates_dT(k_{self.secondary_rates[0].cname()});\n"
            fstring += f"        {dtype} drdT_pg = rate_eval.dscreened_rates_dT(k_{self.secondary_rates[1].cname()});\n"
            fstring += f"        {dtype} drdT_pa = rate_eval.dscreened_rates_dT(k_{self.secondary_reverse[1].cname()});\n"
            fstring += "        drate_dT = drdT_ag + drdT_ap * r_pg * dd + r_ap * drdT_pg * dd - r_ap * r_pg * dd * dd * (drdT_pg + drdT_pa);\n"
            fstring += "    }\n"
        else:

            # first we need to get all of the rates that make this up
            fstring += f"    {dtype} r_ga = rate_eval.screened_rates(k_{self.primary_reverse.cname()});\n"
            fstring += f"    {dtype} r_pa = rate_eval.screened_rates(k_{self.secondary_reverse[1].cname()});\n"
            fstring += f"    {dtype} r_gp = rate_eval.screened_rates(k_{self.secondary_reverse[0].cname()});\n"
            fstring += f"    {dtype} r_pg = rate_eval.screened_rates(k_{self.secondary_rates[1].cname()});\n"

            # now the approximation
            fstring += f"    {dtype} dd = 1.0_rt / (r_pg + r_pa);\n"
            fstring += "    rate = r_ga + r_gp * r_pa * dd;\n"
            fstring += "    if constexpr (std::is_same_v<T, rate_derivs_t>) {\n"
            fstring += f"        {dtype} drdT_ga = rate_eval.dscreened_rates_dT(k_{self.primary_reverse.cname()});\n"
            fstring += f"        {dtype} drdT_pa = rate_eval.dscreened_rates_dT(k_{self.secondary_reverse[1].cname()});\n"
            fstring += f"        {dtype} drdT_gp = rate_eval.dscreened_rates_dT(k_{self.secondary_reverse[0].cname()});\n"
            fstring += f"        {dtype} drdT_pg = rate_eval.dscreened_rates_dT(k_{self.secondary_rates[1].cname()});\n"
            fstring += "        drate_dT = drdT_ga + drdT_gp * r_pa * dd + r_gp * drdT_pa * dd - r_gp * r_pa * dd * dd * (drdT_pg + drdT_pa);\n"
            fstring += "    }\n"

        if not leave_open:
            fstring += "}\n\n"

        return fstring
