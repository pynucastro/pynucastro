"""Classes and method for describing effective rates that can be
limited by a beta waiting point in the sequence.

"""

import copy

import numpy as np

from pynucastro.nucdata import Nucleus
from pynucastro.rates.rate import Rate
from pynucastro.rates.reaclib_rate import ReacLibRate
from pynucastro.rates.starlib_rate import StarLibRate
from pynucastro.rates.temperature_tabular_rate import TemperatureTabularRate


class BetaLimitedRate(Rate):
    """A beta-limited rate represents a sequence that at low
    temperature follows the behavior of the underlying rate but at
    high temperature, one or more beta-decays limit the sequence,
    effectively capping the rate.

    An example is the CNO cycle sequence: N14(p,γ)O15(β+)N15.  At low
    temperatures, N14(p,γ)O15 is the limiting rate, but at higher
    temperatures, O15(β+)N15 becomes the waiting point.  So we want to
    cap this rate at the rate for that decay.  Additionally, for hot-CNO,
    we'll have to wait on O14(β+)N14, so we accept a list of weak rates
    and compute the total waiting time from the sequence as the cap.

    Note: a ``BetaLimitedRate`` does not change the products or
    stoichiometry.  It is essentially just a wrapper around a rate
    that provides a cap to the reaction rate.  The idea is that
    a ``BetaLimitedRate`` can be composed with a ``ModifiedRate``
    or ``BranchedRate`` as needed to do more complex operations.

    Parameters
    ----------
    underlying_rate : Rate
        the underlying rate we are evaluating numerically to get the
        number of reactions / sec if there were no beta-limiting.
    beta_limiting_rates : list(Rate)
    description : str
        a description of the rate sequence we are approximating.  This
        will be added as a comment to code outputs.

    """

    def __init__(self, underlying_rate, beta_limiting_rates, *,
                 limiter_nucleus=None, description=None):

        self.underlying_rate = underlying_rate

        if not isinstance(beta_limiting_rates, (list, tuple)):
            beta_limiting_rates = list(beta_limiting_rates)
        self.beta_limiting_rates = beta_limiting_rates

        self.limiter_nucleus = Nucleus.cast(limiter_nucleus)
        self.description = description

        # we expect an unmodified rate as the underlying rate
        assert isinstance(underlying_rate,
                          (ReacLibRate, StarLibRate, TemperatureTabularRate))

        # we also assume that the underlying rate is a two-body
        # reaction, whose flux will take the form F =
        # ρY(ξ)Y(X)λ_{X(ξ,γ)}
        assert underlying_rate.dens_exp == 1

        # also make sure that limiter nucleus is a reactant of the
        # underlying rate
        assert self.limiter_nucleus in underlying_rate.reactants

        # make sure that the beta-limiting rates are weak rates
        assert all(r.weak_type in ["beta_pos", "beta_neg"]
                   for r in beta_limiting_rates)

        # the reactants and products are the same as the underlying
        # rate
        reactants = self.underlying_rate.reactants
        products = self.underlying_rate.products

        # we will need composition to be able to compare rates
        self.rate_eval_needs_rho = True
        self.rate_eval_needs_comp = True

        super().__init__(reactants=reactants, products=products,
                         label="betalimited")

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

        # copy the original rate
        new.underlying_rate = copy.copy(self.underlying_rate)

        # and the beta rates
        new.beta_limiting_rates = [copy.copy(r)
                                   for r in self.beta_limiting_rates]

        return new

    def get_child_rates(self):
        """Return a list of all of the rates that are used in this
        approximation.

        Returns
        -------
        list(Rate)

        """
        return [self.underlying_rate, *self.beta_limiting_rates]

    def log_eval(self, T, *, rho=None, comp=None,
                 screen_func=None):
        """Evaluate the natural log of reaction rate for the
        beta-limited rate.

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

        return np.log(self.eval(T, rho=rho, comp=comp,
                                screen_func=screen_func))

    def eval(self, T, *, rho=None, comp=None,
             screen_func=None):
        """Evaluate the beta-limited rate.

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

        # evaluate the underlying rate -- this will just be
        # N_A <σv>
        r0 = self.underlying_rate.eval(T, rho=rho, comp=comp,
                                       screen_func=screen_func)

        # evaluate the beta-limiting effective rate
        lambdas = [rate.eval(T, rho=rho, comp=comp)
                   for rate in self.beta_limiting_rates]

        # sequential beta waiting times add.
        lambda_beta_tot = 1.0 / sum(1.0 / lam for lam in lambdas)

        # the other reactant in the underlying rate reaction
        Y_limiter = comp.get_molar()[self.limiter_nucleus]

        # the core idea here is that the underlying rate flux is evaluated
        # as F = ρY(ξ)Y(X)λ_{X(ξ,γ)} and the β-limiting rate is evaluated
        # simply as F = Y(X)λ_β, and for the final rate, we want
        # F = ρY(ξ)Y(X) min{λ_{X(ξ,γ)}, λ_β/(ρY(ξ))}

        return min(r0, lambda_beta_tot / (rho * Y_limiter))
