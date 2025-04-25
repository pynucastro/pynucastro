from pynucastro.nucdata import Nucleus
from pynucastro.rates.rate import Rate


class ModifiedRate(Rate):
    """A modified rate takes an original rate and changes some
    properties of it.  The evaluation of the original rate will still
    be used for the actual rate, but the modified rate can have a
    different products (and therefore Q value) or stoichiometric
    coefficients

    Parameters
    ----------
    original_rate : Rate
        the underlying rate we are evaluating numerically to
        get the number of reactions / sec (with suitable volume
        scalings)
    stoichiometry : dict(Nucleus)
        a custon set of coefficients to be used in the
        evolution equations dY(Nucleus)/dt.  If this is not
        set, then simply the count of each nucleus in the
        list of reactants and products will be used.
    new_products : list(Nucleus)
        a list of nuclei that should be used as the product
        of the modified rate, instead of the products from the
        original rate.
    """

    def __init__(self, original_rate, *,
                 stoichiometry=None,
                 new_products=None):

        self.original_rate = original_rate

        reactants = original_rate.reactants
        if new_products:
            products = new_products
        else:
            products = original_rate.products

        super().__init__(reactants=reactants, products=products,
                         label="modified",
                         stoichiometry=stoichiometry)

        self.chapter = "m"

        # update the Q value
        if new_products:
            self._set_q()

    def eval(self, T, *, rho=None, comp=None):
        """evaluate the approximate rate"""

        if self.approx_type == "ap_pg":
            if not self.is_reverse:  # pylint: disable=no-else-return
                # the approximate forward rate is r_ag + r_ap r_pg / (r_pg + r_pa)
                r_ag = self.rates["A(a,g)B"].eval(T)
                r_ap = self.rates["A(a,p)X"].eval(T)
                r_pg = self.rates["X(p,g)B"].eval(T)

                r_pa = self.rates["X(p,a)A"].eval(T)

                return r_ag + r_ap * r_pg / (r_pg + r_pa)

            else:
                # the approximate reverse rate is r_ga + r_pa r_gp / (r_pg + r_pa)

                r_ga = self.rates["B(g,a)A"].eval(T)
                r_gp = self.rates["B(g,p)X"].eval(T)
                r_pa = self.rates["X(p,a)A"].eval(T)

                r_pg = self.rates["X(p,g)B"].eval(T)

                return r_ga + r_pa * r_gp / (r_pg + r_pa)

        elif self.approx_type == "nn_g":

            # we are approximating A(n,g)X(n,g)B

            Yn = comp.get_molar()[Nucleus("n")]

            if not self.is_reverse:  # pylint: disable=no-else-return
                # the forward rate
                A_ng_X = self.rates["A(n,g)X"].eval(T)
                X_ng_B = self.rates["X(n,g)B"].eval(T)

                X_gn_A = self.rates["X(g,n)A"].eval(T)

                return A_ng_X * X_ng_B / (rho * Yn * X_ng_B + X_gn_A)

            else:
                # the reverse rate
                B_gn_X = self.rates["B(g,n)X"].eval(T)
                X_gn_A = self.rates["X(g,n)A"].eval(T)

                X_ng_B = self.rates["X(n,g)B"].eval(T)

                return B_gn_X * X_gn_A / (rho * Yn * X_ng_B + X_gn_A)

        raise NotImplementedError(f"approximation type {self.approx_type} not supported")

    def function_string_py(self):
        """
        Return a string containing python function that computes the
        approximate rate
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

    def function_string_cxx(self, dtype="double", specifiers="inline", leave_open=False, extra_args=()):
        """
        Return a string containing C++ function that computes the
        approximate rate
        """

        if dtype == "amrex::Real":
            array_type = "amrex::Array1D"
        else:
            array_type = "Array1D"

        if self.approx_type == "ap_pg":

            args = ["const T& rate_eval", f"{dtype}& rate", f"{dtype}& drate_dT", *extra_args]
            fstring = ""
            fstring = "template <typename T>\n"
            fstring += f"{specifiers}\n"
            fstring += f"void rate_{self.cname()}({', '.join(args)}) {{\n\n"

            if not self.is_reverse:

                # first we need to get all of the rates that make this up
                fstring += f"    {dtype} r_ag = rate_eval.screened_rates(k_{self.rates['A(a,g)B'].cname()});\n"
                fstring += f"    {dtype} r_ap = rate_eval.screened_rates(k_{self.rates['A(a,p)X'].cname()});\n"
                fstring += f"    {dtype} r_pg = rate_eval.screened_rates(k_{self.rates['X(p,g)B'].cname()});\n"
                fstring += f"    {dtype} r_pa = rate_eval.screened_rates(k_{self.rates['X(p,a)A'].cname()});\n"

                # now the approximation
                fstring += f"    {dtype} dd = 1.0_rt / (r_pg + r_pa);\n"
                fstring += "    rate = r_ag + r_ap * r_pg * dd;\n"
                fstring += "    if constexpr (std::is_same_v<T, rate_derivs_t>) {\n"
                fstring += f"        {dtype} drdT_ag = rate_eval.dscreened_rates_dT(k_{self.rates['A(a,g)B'].cname()});\n"
                fstring += f"        {dtype} drdT_ap = rate_eval.dscreened_rates_dT(k_{self.rates['A(a,p)X'].cname()});\n"
                fstring += f"        {dtype} drdT_pg = rate_eval.dscreened_rates_dT(k_{self.rates['X(p,g)B'].cname()});\n"
                fstring += f"        {dtype} drdT_pa = rate_eval.dscreened_rates_dT(k_{self.rates['X(p,a)A'].cname()});\n"
                fstring += "        drate_dT = drdT_ag + drdT_ap * r_pg * dd + r_ap * drdT_pg * dd - r_ap * r_pg * dd * dd * (drdT_pg + drdT_pa);\n"
                fstring += "    }\n"
            else:

                # first we need to get all of the rates that make this up
                fstring += f"    {dtype} r_ga = rate_eval.screened_rates(k_{self.rates['B(g,a)A'].cname()});\n"
                fstring += f"    {dtype} r_pa = rate_eval.screened_rates(k_{self.rates['X(p,a)A'].cname()});\n"
                fstring += f"    {dtype} r_gp = rate_eval.screened_rates(k_{self.rates['B(g,p)X'].cname()});\n"
                fstring += f"    {dtype} r_pg = rate_eval.screened_rates(k_{self.rates['X(p,g)B'].cname()});\n"

                # now the approximation
                fstring += f"    {dtype} dd = 1.0_rt / (r_pg + r_pa);\n"
                fstring += "    rate = r_ga + r_gp * r_pa * dd;\n"
                fstring += "    if constexpr (std::is_same_v<T, rate_derivs_t>) {\n"
                fstring += f"        {dtype} drdT_ga = rate_eval.dscreened_rates_dT(k_{self.rates['B(g,a)A'].cname()});\n"
                fstring += f"        {dtype} drdT_pa = rate_eval.dscreened_rates_dT(k_{self.rates['X(p,a)A'].cname()});\n"
                fstring += f"        {dtype} drdT_gp = rate_eval.dscreened_rates_dT(k_{self.rates['B(g,p)X'].cname()});\n"
                fstring += f"        {dtype} drdT_pg = rate_eval.dscreened_rates_dT(k_{self.rates['X(p,g)B'].cname()});\n"
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
            fstring += f"void rate_{self.cname()}({', '.join(args)}) {{\n\n"
            fstring += f"    {dtype} Yn = Y(N);\n"

            if not self.is_reverse:

                # first we need to get all of the rates that make this up
                fstring += f"    {dtype} r1_ng = rate_eval.screened_rates(k_{self.rates['A(n,g)X'].cname()});\n"
                fstring += f"    {dtype} r2_ng = rate_eval.screened_rates(k_{self.rates['X(n,g)B'].cname()});\n"
                fstring += f"    {dtype} r1_gn = rate_eval.screened_rates(k_{self.rates['X(g,n)A'].cname()});\n"

                # now the approximation
                fstring += f"    {dtype} dd = 1.0_rt / (rho * Yn * r2_ng + r1_gn);\n"
                fstring += "    rate = r1_ng * r2_ng * dd;\n"
                fstring += "    if constexpr (std::is_same_v<T, rate_derivs_t>) {\n"
                fstring += f"        {dtype} dr1dT_ng = rate_eval.dscreened_rates_dT(k_{self.rates['A(n,g)X'].cname()});\n"
                fstring += f"        {dtype} dr2dT_ng = rate_eval.dscreened_rates_dT(k_{self.rates['X(n,g)B'].cname()});\n"
                fstring += f"        {dtype} dr1dT_gn = rate_eval.dscreened_rates_dT(k_{self.rates['X(g,n)A'].cname()});\n"
                fstring += "        drate_dT = dr1dT_ng * r2_ng * dd + r1_ng * dr2dT_ng * dd - r1_ng * r2_ng * dd * dd * (rho * Yn * dr2dT_ng + dr1dT_gn);\n"
                fstring += "    }\n"
            else:

                # first we need to get all of the rates that make this up
                fstring += f"    {dtype} r1_gn = rate_eval.screened_rates(k_{self.rates['X(g,n)A'].cname()});\n"
                fstring += f"    {dtype} r2_gn = rate_eval.screened_rates(k_{self.rates['B(g,n)X'].cname()});\n"
                fstring += f"    {dtype} r2_ng = rate_eval.screened_rates(k_{self.rates['X(n,g)B'].cname()});\n"

                # now the approximation
                fstring += f"    {dtype} dd = 1.0_rt / (rho * Yn * r2_ng + r1_gn);\n"
                fstring += "    rate = r1_gn * r2_gn * dd;\n"
                fstring += "    if constexpr (std::is_same_v<T, rate_derivs_t>) {\n"
                fstring += f"        {dtype} dr1dT_gn = rate_eval.dscreened_rates_dT(k_{self.rates['X(g,n)A'].cname()});\n"
                fstring += f"        {dtype} dr2dT_gn = rate_eval.dscreened_rates_dT(k_{self.rates['B(g,n)X'].cname()});\n"
                fstring += f"        {dtype} dr2dT_ng = rate_eval.dscreened_rates_dT(k_{self.rates['X(n,g)B'].cname()});\n"
                fstring += "        drate_dT = dr1dT_gn * r2_gn * dd + r1_gn * dr2dT_gn * dd - r1_gn * r2_gn * dd * dd * (rho * Yn * dr2dT_ng + dr1dT_gn);\n"
                fstring += "    }\n"

            if not leave_open:
                fstring += "}\n\n"

            return fstring

        raise NotImplementedError("don't know how to work with this approximation")
