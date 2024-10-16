from pynucastro.nucdata import Nucleus
from pynucastro.rates.rate import Rate


class ApproximateRate(Rate):

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
                                 products=[self.primary_product], label="approx")
            else:
                super().__init__(reactants=[self.primary_product],
                                 products=[self.primary_reactant, Nucleus("he4")], label="approx")

            self.chapter = "a"

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
