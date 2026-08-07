"""Classes and methods for describing rates where one or more
properties have been modified from the original source.

"""

import copy

import numpy as np

from pynucastro.rates.rate import Rate, ThermoState
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

       N14(p,γ)O15(,e⁺ν)N15(p,α)C12
       N14(p,γ)O15(,e⁺ν)N15(p,γ)O16

    These differ only in the last rate.  We would set the
    underlying_rate to be N14(p,γ)O15, the primary_branch to be
    N15(p,α)C12 and the seconary branch to be N15(p,γ))O16.  It would
    then compute the branching ratio:

       f = λ_{N15(p,α)C12} / (λ_{N15(p,α)C12} + λ_{N15(p,γ)O16)

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
    description : str
        a description of the rate sequence we are approximating.  This
        will be added as a comment to code outputs.

    """

    def __init__(self, underlying_rate, *,
                 primary_branch=None,
                 other_branch=None,
                 description=None):

        self.underlying_rate = underlying_rate
        self.primary_branch = primary_branch
        self.other_branch = other_branch

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
                         label="branched")

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
        new.original_rate = copy.copy(self.original_rate)

        return new

    def log_eval(self, T, *, rho=None, comp=None,
                 screen_func=None):
        """Evaluate natural log of reaction rates for the modified rate.
        This simply calls the evaluation of the underlying original rate.

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
        numpy.ndarray

        """

        # Evaluate original rate without screening
        # The modified rate can have a different set of reactants for screening
        log_rate = self.original_rate.log_eval(T, rho=rho, comp=comp, screen_func=None)

        # Apply screening correction
        log_scor = 0.0
        if screen_func is not None:
            if rho is None or comp is None:
                raise ValueError("rho (density) and comp (Composition) needs to be defined when applying electron screening.")
            state = ThermoState(rho=rho, T=T, comp=comp)
            log_scor = self.evaluate_screening(state, screen_func=screen_func)

        # To consider general cases, convert to 1D array
        log_rate = np.atleast_1d(log_rate)

        # Apply screening
        log_rate += log_scor

        return log_rate

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
        fstring += f"def {self.fname}(rate_eval, tf, log_scor=0.0):\n"
        fstring += f"    # {self.rid}\n"
        fstring += f"    {self.original_rate.fname}(rate_eval, tf, log_scor=log_scor)\n"
        fstring += f"    rate_eval.{self.fname} = rate_eval.{self.original_rate.fname}\n\n"
        return fstring

    def function_string_cxx(self, dtype="double", specifiers="inline",
                            leave_open=False, extra_args=()):
        """Return a string containing the C++ function that computes
        the rate.  For a ModifiedRate, this simply calls the
        corresponding function for the underlying original rate.

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

        args = ["const tf_t& tfactors",
                f"const {dtype} log_scor", f"const {dtype} dlog_scor_dT",
                f"{dtype}& rate", f"{dtype}& drate_dT", *extra_args]
        fstring = ""
        fstring = "template <int do_T_derivatives>\n"
        fstring += f"{specifiers}\n"
        fstring += f"void rate_{self.fname}({', '.join(args)}) {{\n\n"

        # first we need to get all of the rates that make this up
        fstring += f"    // {self.rid} (calls the underlying rate)\n\n"
        fstring += f"    rate_{self.original_rate.fname}<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);\n"

        if not leave_open:
            fstring += "}\n\n"

        return fstring
