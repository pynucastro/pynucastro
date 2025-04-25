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
        a custom set of coefficients to be used in the
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

        self.original_rate.eval(T, rho=rho, comp=comp)

    def function_string_py(self):
        """Return a string containing the python function that
        computes the rate -- in this case it is the underlying
        original rate.

        Returns
        -------
        str

        """

        return self.original_rate.function_string_py()

    def function_string_cxx(self, dtype="double", specifiers="inline",
                            leave_open=False, extra_args=()):
        """
        Return a string containing C++ function that computes the
        approximate rate
        """

        return self.function_string_cxx(dtype=dtype, specifiers=specifiers,
                                        leave_open=leave_open, extra_args=extra_args)
