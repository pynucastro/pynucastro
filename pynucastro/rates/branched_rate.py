"""Classes and methods for describing rate sequences that
have branching endpoints.

"""

import copy

import numpy as np

from pynucastro.rates.rate import Rate
from pynucastro.rates.reaclib_rate import ReacLibRate
from pynucastro.rates.starlib_rate import StarLibRate
from pynucastro.rates.temperature_tabular_rate import TemperatureTabularRate


class BranchedRate(Rate):
    """A branched rate represents a sequence that can have different
    endpoints depending on branching.  It takes an underlying_rate
    which will be used to evaluate the rate, and then takes a
    primary_branch and other_branches rate(s) that are used to
    normalize the rate.  The products of the rate are set to be
    the products of the primary_branch.

    An example application would be the sequences:

    * N14(p,γ)O15(,e⁺ν)N15(p,α)C12
    * N14(p,γ)O15(,e⁺ν)N15(p,γ)O16

    These differ only in the last rate.  We would set the
    underlying_rate to be N14(p,γ)O15, the primary_branch to be
    N15(p,α)C12 and the secondary branch to be N15(p,γ))O16.  It would
    then compute the branching ratio:

    f = λ_{N15(p,α)C12} / (λ_{N15(p,α)C12} + λ_{N15(p,γ)O16))

    and the final rate evaluation would be

    λ = f λ_{N14(p,γ)O15}

    Parameters
    ----------
    underlying_rate : Rate
        the underlying rate we are evaluating numerically to get the
        number of reactions / sec (with suitable volume scalings),
        reduced by the branching fraction
    primary_branch : Rate
        the branch we want this sequence to use
    other_branch : Rate
        an alternate branch used in normalization
    stoichiometry : dict(Nucleus)
        a custom set of coefficients to be used in the evolution
        equations dY(Nucleus)/dt.  If this is not set, then simply the
        count of each nucleus in the list of reactants and products
        will be used.
    description : str
        a description of the rate sequence we are approximating.  This
        will be added as a comment to code outputs.

    """

    def __init__(self, underlying_rate, *,
                 primary_branch=None,
                 other_branch=None,
                 stoichiometry=None,
                 description=None):

        self.underlying_rate = underlying_rate
        self.primary_branch = primary_branch
        self.other_branch = other_branch
        self.description = description

        # at the moment, this is only tested with ReacLibRate,
        # TemperatureTabularRate, and StarLibRate rates.  It is
        # important in the C++ code generation the we fill branched rates
        # only after the other rates are filled and screened.

        assert isinstance(underlying_rate,
                          (ReacLibRate, StarLibRate, TemperatureTabularRate))

        assert isinstance(primary_branch,
                          (ReacLibRate, StarLibRate, TemperatureTabularRate))

        assert isinstance(other_branch,
                          (ReacLibRate, StarLibRate, TemperatureTabularRate))

        # the reactants come from the underlying rate
        reactants = self.underlying_rate.reactants

        # the products come from the branch we take
        products = self.primary_branch.products

        super().__init__(reactants=reactants, products=products,
                         weak_type=self.underlying_rate.weak_type,
                         stoichiometry=stoichiometry,
                         label="branched")

        # for the moment, we only work if both branches have the same
        # reactants.  If they don't then we need to weight by (rho Y)
        # for each nucleus they don't have in common.  We'll also
        # need to set rate_eval_needs_rho and rate_eval_needs_comp
        assert self.primary_branch.reactants == self.other_branch.reactants

        self._set_print_representation()

    def __copy__(self):
        """Make a copy of the rate via copy.copy().  This is mostly
        shallow except for a few attributes to address some mutability
        issues

        """

        cls = type(self)
        new = cls.__new__(cls)

        # shallow copy everything
        new.__dict__ = self.__dict__.copy()

        # override some shallow copies
        new.reactants = list(self.reactants)
        new.products = list(self.products)
        if self.stoichiometry:
            new.stoichiometry = dict(self.stoichiometry)

        # copy the original rate
        new.underlying_rate = copy.copy(self.underlying_rate)

        return new

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

        return np.log(self.eval(T, rho=rho, comp=comp, screen_func=screen_func))

    def eval(self, T, *, rho=None, comp=None,
             screen_func=None):
        """Evaluate the branched rate.

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

        # evaluate the underlying rate
        r0 = self.underlying_rate.eval(T, rho=rho, comp=comp,
                                       screen_func=screen_func)

        # now evaluate the branches
        r_prim_br = self.primary_branch.eval(T, rho=rho, comp=comp,
                                           screen_func=screen_func)
        r_other_br = self.other_branch.eval(T, rho=rho, comp=comp,
                                          screen_func=screen_func)

        # compute the branching factor
        f = r_prim_br / (r_prim_br + r_other_br)

        return f * r0

    def get_child_rates(self):
        """Return a list of all of the rates that are used in this
        approximation.

        Returns
        -------
        list(Rate)

        """
        return [self.underlying_rate, self.primary_branch, self.other_branch]

    def function_string_py(self):
        """Return a string containing the python function that
        computes the rate -- in this case it is the underlying rate
        modified by the branching ratio

        Returns
        -------
        str

        """

        fstring = ""
        fstring += "@numba.njit()\n"
        fstring += f"def {self.fname}(rate_eval, tf, log_scor=0.0):\n"
        fstring += f"    # {self.rid}\n"
        if self.description:
            fstring += f"    # represents the sequence: {self.description}\n\n"
        fstring += f"    r0 = rate_eval.{self.underlying_rate.fname}\n"
        fstring += f"    r_prim_br = rate_eval.{self.primary_branch.fname}\n"
        fstring += f"    r_other_br = rate_eval.{self.other_branch.fname}\n\n"
        fstring += "    f = r_prim_br / (r_prim_br + r_other_br)\n"
        fstring += f"    rate_eval.{self.fname} = f * r0\n\n"
        return fstring

    def function_string_cxx(self, dtype="double", specifiers="inline",
                            leave_open=False, extra_args=()):
        """Return a string containing the C++ function that computes
        the rate.  For a BranchedRate, this returns the underlying
        original rate modified by the branching ratio.

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
        extra_args : list, tuple
            A list of strings representing additional arguments that
            should be appended to the argument list when defining the
            function interface.

        Returns
        -------
        str

        """

        args = ["const T& rate_eval", f"{dtype}& rate", f"{dtype}& drate_dT", *extra_args]
        fstring = ""
        fstring = "template <typename T>\n"
        fstring += f"{specifiers}\n"
        fstring += f"void rate_{self.fname}({', '.join(args)}) {{\n\n"

        fstring += f"    // {self.rid} (branched rate)\n\n"
        if self.description:
            fstring += f"    // represents the sequence: {self.description}\n\n"

        fstring += f"    {dtype} r0 = rate_eval.screened_rates(k_{self.underlying_rate.fname});\n"
        fstring += f"    {dtype} r_prim_br = rate_eval.screened_rates(k_{self.primary_branch.fname});\n"
        fstring += f"    {dtype} r_other_br = rate_eval.screened_rates(k_{self.other_branch.fname});\n\n"

        # now do the approximation
        fstring += f"    {dtype} f = r_prim_br / (r_prim_br + r_other_br);\n"
        fstring += "    rate = f * r0;\n\n"

        fstring += "    if constexpr (std::is_same_v<T, rate_derivs_t>) {\n"
        fstring += f"        {dtype} drdT_0 = rate_eval.dscreened_rates_dT(k_{self.underlying_rate.fname});\n"
        fstring += f"        {dtype} drdT_prim_br = rate_eval.dscreened_rates_dT(k_{self.primary_branch.fname});\n"
        fstring += f"        {dtype} drdT_other_br = rate_eval.dscreened_rates_dT(k_{self.other_branch.fname});\n\n"

        fstring += f"        {dtype} dfdT = (drdT_prim_br - f * (drdT_prim_br + drdT_other_br)) / (r_prim_br + r_other_br);\n"
        fstring += "        drate_dT = f * drdT_0 + dfdT * r0;\n"
        fstring += "    }\n"

        if not leave_open:
            fstring += "}\n\n"

        return fstring
