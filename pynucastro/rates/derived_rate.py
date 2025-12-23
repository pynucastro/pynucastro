"""Classes and methods for deriving inverse rates from forward rates
using detailed balance.

"""

import math
import warnings
from collections import Counter

import numpy as np

from pynucastro.constants import constants
from pynucastro.nucdata import Nucleus
from pynucastro.rates.rate import Rate, Tfactors
from pynucastro.rates.reaclib_rate import ReacLibRate, SingleSet
from pynucastro.rates.tabular_rate import TabularRate


class DerivedRate(Rate):
    """A reverse rate computed from a forward rate via detailed
    balance.

    Parameters
    ----------
    source_rate : Rate
        The forward rate that will be used to derive the reverse
    use_pf : bool
        Do we apply the partition function?
    use_unreliable_spins : bool
        Do we use spins that are weakly experimentally supported?

    """

    def __init__(self, source_rate, use_pf=False, use_unreliable_spins=True):

        self.use_pf = use_pf
        self.source_rate = source_rate
        self.use_unreliable_spins = use_unreliable_spins

        if not isinstance(self.source_rate, Rate):
            raise TypeError('The source rate must be a Rate subclass')

        if (isinstance(self.source_rate, TabularRate) or self.source_rate.weak or
            self.source_rate.derived_from_inverse):
            raise ValueError('The source rate is a ReacLib derived from inverse rate or weak or tabular')

        if self.source_rate.stoichiometry is not None:
            warnings.warn(UserWarning('Using rates with stoichiometry will not be compatible with NSE.'))

        if not all(nuc.spin_states for nuc in self.source_rate.reactants):
            raise ValueError(f'One of the reactants spin ground state ({self.source_rate.reactants}), is not defined')

        if not all(nuc.spin_states for nuc in self.source_rate.products):
            raise ValueError(f'One of the products spin ground state ({self.source_rate.products}), is not defined')

        if not use_unreliable_spins:
            if not all(nuc.spin_reliable for nuc in self.source_rate.reactants):
                raise ValueError(f'One of the reactants spin ground state ({self.source_rate.reactants}), is considered unreliable')
            if not all(nuc.spin_reliable for nuc in self.source_rate.products):
                raise ValueError(f'One of the products spin ground state ({self.source_rate.products}), is considered unreliable')

        # Compute temperature-independent prefactor of the equilibrium ratio
        F = 1.0

        F *= math.prod(nucr.spin_states for nucr in self.source_rate.reactants)
        F /= math.prod(nucr.A for nucr in self.source_rate.reactants)
        F *= math.prod(nucr.A_nuc for nucr in self.source_rate.reactants)**2.5

        F /= math.prod(nucp.spin_states for nucp in self.source_rate.products)
        F *= math.prod(nucp.A for nucp in self.source_rate.products)
        F /= math.prod(nucp.A_nuc for nucp in self.source_rate.products)**2.5

        F *= self.counter_factors()[1] / self.counter_factors()[0]

        self.net_stoich = len(self.source_rate.reactants) - len(self.source_rate.products)
        if self.net_stoich != 0:
            F *= constants.m_u_C18**(2.5 * self.net_stoich)
            F *= (constants.k / (2.0 * np.pi * constants.hbar**2))**(1.5 * self.net_stoich)

        self.ratio_factor = F

        super().__init__(reactants=self.source_rate.products,
                         products=self.source_rate.reactants,
                         label="derived", rate_source=self.source_rate.src,
                         stoichiometry=self.source_rate.stoichiometry)

        # If source rate is a reaclib rate, then create a derived reaclib set based
        # on the source rate by absorbing the equilibrium ratio within the
        # reaclib log terms for better numerical stability at lower temp
        self.derived_sets = None
        if isinstance(self.source_rate, ReacLibRate):
            self.derived_sets = []
            for source_set in self.source_rate.sets:
                a_derived = source_set.a.copy()
                a_derived[0] += np.log(self.ratio_factor) + 13.5 * self.net_stoich * np.log(10)
                a_derived[1] += self.Q / (constants.k_MeV * 1.0e9)
                a_derived[6] += 1.5 * self.net_stoich
                self.derived_sets.append(SingleSet(a_derived, source_set.labelprops))

        # explicitly mark it as a rate derived from inverse
        self.derived_from_inverse = True

    def _warn_about_missing_pf_tables(self):
        skip_nuclei = {Nucleus("h1"), Nucleus("n"), Nucleus("he4")}
        for nuc in set(self.source_rate.reactants + self.source_rate.products) - skip_nuclei:
            if not nuc.partition_function:
                warnings.warn(UserWarning(f'{nuc} partition function is not supported by tables: set pf = 1.0 by default'))

    def eval(self, T, *, rho=None, comp=None,
             screen_func=None):
        """Evaluate the derived reverse rate.

        Parameters
        ----------
        T : float
            the temperature to evaluate the rate at
        rho : float
            the density to evaluate screening effects at.
        comp : float
            the composition (of type
            :py:class:`Composition <pynucastro.networks.rate_collection.Composition>`)
            to evaluate screening effects with.
        screen_func : Callable
            one of the screening functions from :py:mod:`pynucastro.screening`
            -- if provided, then the rate will include screening correction.

        Returns
        -------
        float

        """

        # evaluate the source rate without screening

        if self.derived_sets is not None:
            r = 0.0
            tf = Tfactors(T)
            for s in self.derived_sets:
                f = s.f()
                r += f(tf)
        else:
            r = self.source_rate.eval(T=T, rho=rho, comp=comp, screen_func=None)

            # compute the equilibrium ratio
            ratio = self.ratio_factor * T**(1.5 * self.net_stoich)

            # Note Q-value here is for the derived rate.
            ratio *= np.exp(self.Q / (constants.k_MeV * T))

            # apply ratio
            r *= ratio

        z_r = 1.0
        z_p = 1.0
        if self.use_pf:
            self._warn_about_missing_pf_tables()

            for nucr in self.source_rate.reactants:
                if nucr.partition_function is not None:
                    z_r *= nucr.partition_function.eval(T)

            for nucp in self.source_rate.products:
                if nucp.partition_function is not None:
                    z_p *= nucp.partition_function.eval(T)

        # Apply equilibrium ratio
        r *= z_r / z_p

        # Apply screening correction
        scor = 1.0
        if screen_func is not None:
            if rho is None or comp is None:
                raise ValueError("rho (density) and comp (Composition) needs to be defined when applying electron screening.")
            scor = self.evaluate_screening(rho, T, comp, screen_func)

        r *= scor

        return r

    def function_string_py(self):
        """Return a string containing the python function that
        computes the rate.

        Returns
        -------
        str

        """

        fstring = ""
        fstring += "@numba.njit()\n"
        fstring += f"def {self.fname}(rate_eval, tf):\n"
        fstring += f"    # {self.rid}\n\n"

        if self.derived_sets is not None:
            fstring += "    rate = 0.0 \n\n"
            for s in self.derived_sets:
                fstring += f"    # {s.labelprops[0:5]}\n"
                set_string = s.set_string_py(prefix="rate", plus_equal=True)
                for t in set_string.split("\n"):
                    fstring += "    " + t + "\n"
            fstring += "\n"
            fstring += f"    rate_eval.{self.fname} = rate\n"
        else:
            fstring += "    # Evaluate the equilibrium ratio\n"
            fstring += f"    ratio = {self.ratio_factor} * (tf.T9 * 1.0e9)**(1.5 * {self.net_stoich})\n"
            fstring += f"    ratio *= np.exp({self.Q} / (constants.k_MeV * tf.T9 * 1.0e9))\n"
            fstring += f"    rate_eval.{self.fname} = rate_eval.{self.source_rate.fname} * ratio\n"

        if self.use_pf:
            self._warn_about_missing_pf_tables()

            fstring += "\n"
            for nuc in set(self.source_rate.reactants + self.source_rate.products):
                if nuc.partition_function is not None:
                    fstring += f"    # interpolating {nuc} partition function\n"
                    fstring += f"    {nuc}_pf_exponent = np.interp(tf.T9, xp={nuc}_temp_array, fp=np.log10({nuc}_pf_array))\n"
                    fstring += f"    {nuc}_pf = 10.0**{nuc}_pf_exponent\n"
                else:
                    fstring += f"    # setting {nuc} partition function to 1.0 by default, independent of T\n"
                    fstring += f"    {nuc}_pf = 1.0\n"
                fstring += "\n"

            fstring += "    "
            fstring += "z_r = "
            fstring += "*".join([f"{nucr}_pf" for nucr in self.source_rate.reactants])

            fstring += "\n"
            fstring += "    "
            fstring += "z_p = "
            fstring += "*".join([f"{nucp}_pf" for nucp in self.source_rate.products])

            fstring += "\n"
            fstring += f"    rate_eval.{self.fname} *= z_r / z_p\n\n"

        return fstring

    def function_string_cxx(self, dtype="double", specifiers="inline",
                            leave_open=False, extra_args=None):
        """Return a string containing the C++ function that computes
        the derived reverse rate

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

        args = ["const tf_t& tfactors", f"{dtype}& rate", f"{dtype}& drate_dT",
                "[[maybe_unused]] const T& rate_eval",
                "[[maybe_unused]] part_fun::pf_cache_t& pf_cache", *extra_args]

        fstring = ""
        fstring += "template <typename T>\n"
        fstring += f"{specifiers}\n"
        fstring += f"void rate_{self.fname}({', '.join(args)}) {{\n\n"
        fstring += f"    // {self.rid}\n\n"
        fstring += "    rate = 0.0;\n"
        fstring += "    drate_dT = 0.0;\n\n"

        if self.derived_sets is not None:
            fstring += f"    {dtype} ln_set_rate{{0.0}};\n"
            fstring += f"    {dtype} dln_set_rate_dT9{{0.0}};\n"
            fstring += f"    {dtype} set_rate{{0.0}};\n\n"

            for s in self.derived_sets:
                fstring += f"    // {s.labelprops[0:5]}\n"
                set_string = s.set_string_cxx(prefix="ln_set_rate", plus_equal=False, with_exp=False)
                for t in set_string.split("\n"):
                    fstring += "    " + t + "\n"
                fstring += "\n"

                fstring += "    if constexpr (std::is_same_v<T, rate_derivs_t>) {\n"
                dln_set_string_dT9 = s.dln_set_string_dT9_cxx(prefix="dln_set_rate_dT9", plus_equal=False)
                for t in dln_set_string_dT9.split("\n"):
                    fstring += "        " + t + "\n"
                fstring += "    }\n"
                fstring += "\n"

                fstring += "    // avoid underflows by zeroing rates in [0.0, 1.e-100]\n"
                fstring += "    ln_set_rate = std::max(ln_set_rate, -230.0);\n"
                fstring += "    set_rate = std::exp(ln_set_rate);\n"

                fstring += "    rate += set_rate;\n"

                fstring += "    if constexpr (std::is_same_v<T, rate_derivs_t>) {\n"
                fstring += "        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;\n"
                fstring += "    }\n\n"
        else:
            fstring += "    // Evaluate the equilibrium ratio without partition function\n"
            fstring += f"    {dtype} ratio = {self.ratio_factor};\n"
            fstring += f"    {dtype} Q_kBT = {self.Q} / (C::k_MeV * tfactors.T9 * 1.0e9_rt);\n"
            fstring += "    ratio *= std::exp(Q_kBT);\n"
            if self.net_stoich != 0:
                fstring += f"    ratio *= std::sqrt(amrex::Math::powi<{3 * self.net_stoich}>(tfactors.T9 * 1.0e9_rt));\n\n"

            fstring += "    // Apply the ratio without partition function\n"
            fstring += "    // Note that screening is not yet applied to the inverse rate\n"
            fstring += f"    rate = rate_eval.screened_rates(k_{self.source_rate.fname});\n"

            fstring += "    if constexpr (std::is_same_v<T, rate_derivs_t>) {\n"
            fstring += f"        {dtype} dratio_dT = ratio * tfactors.T9i * 1.0e-9_rt * ({1.5 * self.net_stoich} - Q_kBT);\n"
            fstring += f"        drate_dT = rate_eval.dscreened_rates_dT(k_{self.source_rate.fname});\n"
            fstring += "        drate_dT = drate_dT * ratio + rate * dratio_dT;\n"
            fstring += "    }\n"
            fstring += "    rate *= ratio;\n\n"

        # Right now we have rate and drate_dT without the partition
        # function now the partition function corrections

        if self.use_pf:
            self._warn_about_missing_pf_tables()

            fstring += "    // Now apply partition function effects\n"
            for nuc in set(self.source_rate.reactants + self.source_rate.products):
                fstring += f"    {dtype} {nuc}_pf, d{nuc}_pf_dT;\n"

                if nuc.partition_function:
                    fstring += f"    // interpolating {nuc} partition function\n"
                    fstring += f"    get_partition_function_cached({nuc.cindex()}, tfactors, pf_cache, {nuc}_pf, d{nuc}_pf_dT);\n"
                else:
                    fstring += f"    // setting {nuc} partition function to 1.0 by default, independent of T\n"
                    fstring += f"    {nuc}_pf = 1.0_rt;\n"
                    fstring += f"    d{nuc}_pf_dT = 0.0_rt;\n"
                fstring += "\n"

            fstring += f"    {dtype} z_r = "
            fstring += " * ".join([f"{nucr}_pf" for nucr in self.source_rate.reactants])
            fstring += ";\n"

            fstring += f"    {dtype} z_p = "
            fstring += " * ".join([f"{nucp}_pf" for nucp in self.source_rate.products])
            fstring += ";\n\n"

            # now the derivatives, via chain rule

            fstring += "    if constexpr (std::is_same_v<T, rate_derivs_t>) {\n"

            chain_terms = []
            for n in self.source_rate.reactants:
                chain_terms.append(" * ".join([f"{nucr}_pf" for nucr in self.source_rate.reactants if nucr != n] + [f"d{n}_pf_dT"]))
            fstring += f"        {dtype} dz_r_dT = "
            fstring += " + ".join(chain_terms)
            fstring += ";\n"

            chain_terms = []
            for n in self.source_rate.products:
                chain_terms.append(" * ".join([f"{nucp}_pf" for nucp in self.source_rate.products if nucp != n] + [f"d{n}_pf_dT"]))

            fstring += f"        {dtype} dz_p_dT = "
            fstring += " + ".join(chain_terms)
            fstring += ";\n\n"

            fstring += f"        {dtype} dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);\n\n"

            # final terms

            fstring += "        drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);\n"
            fstring += "    }\n"
            fstring += "    rate *= z_r / z_p;\n\n"

        if not leave_open:
            fstring += "}\n\n"

        return fstring

    def counter_factors(self):
        """Compute the multiplicity factor, nucr! = nucr_1! * ... *
        nucr_r!, for each repeated nucr reactant and nucp! = nucp_1! *
        ... * nucp_p! for each nucp product in a ordered pair (nucr!,
        nucp!). The factors nucr! and nucp! avoid overcounting when
        more than one nuclei is involve in the reaction.  If there is
        no multiplicity, then the factor is 1.0.

        Returns
        -------
        tuple(float, float)

        """

        react_counts = Counter(self.source_rate.reactants)
        prod_counts = Counter(self.source_rate.products)

        reactant_factor = 1.0
        for nuc in set(self.source_rate.reactants):
            reactant_factor *= math.factorial(react_counts[nuc])

        product_factor = 1.0
        for nuc in set(self.source_rate.products):
            product_factor *= math.factorial(prod_counts[nuc])

        return (reactant_factor, product_factor)
