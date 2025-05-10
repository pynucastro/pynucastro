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
        the underlying rate we are evaluating numerically to get the
        number of reactions / sec (with suitable volume scalings)
    stoichiometry : dict(Nucleus)
        a custom set of coefficients to be used in the evolution
        equations dY(Nucleus)/dt.  If this is not set, then simply the
        count of each nucleus in the list of reactants and products
        will be used.
    new_reactants : list(Nucleus)
        a list of nuclei that should be used as the reactants of the
        modified rate, instead of the reactants from the original rate.
    new_products : list(Nucleus)
        a list of nuclei that should be used as the products of the
        modified rate, instead of the products from the original rate.
    update_screening : bool
        do we reset the screening pairs for this rate to reflect any
        new products or stoichiometry? or do we still screen based on
        the underlying rate?

    """

    def __init__(self, original_rate, *,
                 stoichiometry=None,
                 new_reactants=None, new_products=None,
                 update_screening=False):

        self.original_rate = original_rate
        self.update_screening = update_screening

        if new_reactants:
            reactants = new_reactants
        else:
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
        if new_products or new_reactants:
            self._set_q()

    def _set_screening(self):
        """Determine if this rate is eligible for screening and the
        nuclei to use.  In this case, we either use the original rate
        or the modified rate, depending on the value of
        update_screening.

        """
        # Tells if this rate is eligible for screening, and if it is
        # then Rate.ion_screen is a 2-element (3 for 3-alpha) list of
        # Nucleus objects for screening; otherwise it is set to none
        self.ion_screen = []
        if self.update_screening:
            _reac = self.reactants
        else:
            _reac = self.original_rate.reactants
        nucz = [q for q in _reac if q.Z != 0]
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
            if self.update_screening:
                _prod = self.products
            else:
                _prod = self.original_rate.products
            nucz = [q for q in _prod if q.Z != 0]
            if len(nucz) > 1:
                nucz.sort(key=lambda x: x.Z)
                self.symmetric_screen = []
                self.symmetric_screen.append(nucz[0])
                self.symmetric_screen.append(nucz[1])
                if len(nucz) == 3:
                    self.symmetric_screen.append(nucz[2])
        else:
            self.symmetric_screen = self.ion_screen

    def eval(self, T, *, rho=None, comp=None):
        """Evaluate the modified rate.  This simply calls the
        evaluation of the underlying original rate.

        Parameters
        ----------
        T : float
            the temperature to evaluate the rate at
        rho : float
            the density to evaluate the rate at
        comp : Composition
            the composition to evaluate the rate with
        Returns
        -------
        float

        """

        return self.original_rate.eval(T, rho=rho, comp=comp)

    def function_string_py(self):
        """Return a string containing the python function that
        computes the rate -- in this case it is the underlying
        original rate.

        Returns
        -------
        str

        """

        fstring = ""
        fstring += "@numba.njit()\n"
        fstring += f"def {self.fname}(rate_eval, tf):\n"
        fstring += f"    # {self.rid}\n"
        fstring += f"    {self.original_rate.fname}(rate_eval, tf)\n"
        fstring += f"    rate_eval.{self.fname} = rate_eval.{self.original_rate.fname}\n\n"
        return fstring

    def function_string_cxx(self, dtype="double", specifiers="inline",
                            leave_open=False, extra_args=()):
        """
        Return a string containing C++ function that computes the
        approximate rate
        """

        args = ["const tf_t& tfactors", f"{dtype}& rate", f"{dtype}& drate_dT", *extra_args]
        fstring = ""
        fstring = "template <int do_T_derivatives>\n"
        fstring += f"{specifiers}\n"
        fstring += f"void rate_{self.cname()}({', '.join(args)}) {{\n\n"

        # first we need to get all of the rates that make this up
        fstring += f"    // {self.rid} (calls the underlying rate)\n\n"
        fstring += f"    rate_{self.original_rate.cname()}<do_T_derivatives>(tfactors, rate, drate_dT);\n"

        if not leave_open:
            fstring += "}\n\n"

        return fstring
