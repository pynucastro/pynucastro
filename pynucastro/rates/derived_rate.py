import math
import warnings
from collections import Counter

import numpy as np

from pynucastro.constants import constants
from pynucastro.nucdata import Nucleus
from pynucastro.rates.rate import Rate
from pynucastro.rates.reaclib_rate import ReacLibRate, SingleSet
from pynucastro.rates.tabular_rate import TabularRate


class DerivedRate(ReacLibRate):
    """A reverse rate computed from a forward rate via detailed
    balance.

    Parameters
    ----------
    rate : Rate
        The forward rate that will be used to derive the reverse
    compute_Q : bool
        Do we recompute the Q-value of the rate from the masses?
    use_pf : bool
        Do we apply the partition function?

    """

    def __init__(self, rate, compute_Q=False, use_pf=False):

        self.use_pf = use_pf
        self.rate = rate
        self.compute_Q = compute_Q

        if not isinstance(rate, Rate):
            raise TypeError('rate must be a Rate subclass')

        if (isinstance(rate, TabularRate) or self.rate.weak or
            self.rate.reverse):
            raise ValueError('The rate is reverse or weak or tabular')

        if not all(nuc.spin_states for nuc in self.rate.reactants):
            raise ValueError('One of the reactants spin ground state, is not defined')

        if not all(nuc.spin_states for nuc in self.rate.products):
            raise ValueError('One of the products spin ground state, is not defined')

        derived_sets = []
        for ssets in self.rate.sets:
            a = ssets.a
            prefactor = 0.0
            Q = 0.0
            prefactor += -np.log(constants.N_A) * (len(self.rate.reactants) -
                                                   len(self.rate.products))

            for nucr in self.rate.reactants:
                prefactor += 1.5*np.log(nucr.A) + np.log(nucr.spin_states)
                Q += nucr.A_nuc
            for nucp in self.rate.products:
                prefactor += -1.5*np.log(nucp.A) - np.log(nucp.spin_states)
                Q -= nucp.A_nuc

            if self.compute_Q:
                Q = Q * constants.m_u_MeV_C18
            else:
                Q = self.rate.Q

            prefactor += np.log(self.counter_factors()[1]) - np.log(self.counter_factors()[0])

            if len(self.rate.reactants) == len(self.rate.products):
                prefactor += 0.0
            else:
                F = (constants.m_u_C18 * constants.k * 1.0e9 /
                     (2.0*np.pi*constants.hbar**2))**(1.5*(len(self.rate.reactants) -
                                                           len(self.rate.products)))
                prefactor += np.log(F)

            a_rev = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            a_rev[0] = prefactor + a[0]
            a_rev[1] = a[1] - Q / (1.0e9 * constants.k_MeV)
            a_rev[2] = a[2]
            a_rev[3] = a[3]
            a_rev[4] = a[4]
            a_rev[5] = a[5]
            a_rev[6] = a[6] + 1.5*(len(self.rate.reactants) -
                                   len(self.rate.products))
            sset_d = SingleSet(a=a_rev, labelprops=rate.labelprops)
            derived_sets.append(sset_d)

        super().__init__(rfile=self.rate.rfile, chapter=self.rate.chapter,
                         original_source=self.rate.original_source,
                         reactants=self.rate.products,
                         products=self.rate.reactants,
                         sets=derived_sets, labelprops="derived", Q=-Q)

        # explicitly mark it as reverse
        self.reverse = True

    def _warn_about_missing_pf_tables(self):
        skip_nuclei = {Nucleus("h1"), Nucleus("n"), Nucleus("he4")}
        for nuc in set(self.rate.reactants + self.rate.products) - skip_nuclei:
            if not nuc.partition_function:
                warnings.warn(UserWarning(f'{nuc} partition function is not supported by tables: set pf = 1.0 by default'))

    def eval(self, T, *, rho=None, comp=None):
        """Evaluate the derived reverse rate.

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

        r = super().eval(T=T, rho=rho, comp=comp)
        z_r = 1.0
        z_p = 1.0
        if self.use_pf:
            self._warn_about_missing_pf_tables()

            for nucr in self.rate.reactants:
                if not nucr.partition_function:
                    continue
                    #nucr.partition_function = lambda T: 1.0
                z_r *= nucr.partition_function.eval(T)

            for nucp in self.rate.products:
                if not nucp.partition_function:
                    continue
                    #nucp.partition_function = lambda T: 1.0
                z_p *= nucp.partition_function.eval(T)

            return r*z_r/z_p
        return r

    def function_string_py(self):
        """Return a string containing python function that computes
        the rate.

        """

        self._warn_about_missing_pf_tables()

        fstring = super().function_string_py()

        if self.use_pf:

            fstring += "\n"
            for nuc in set(self.rate.reactants + self.rate.products):
                if nuc.partition_function:
                    fstring += f"    # interpolating {nuc} partition function\n"
                    fstring += f"    {nuc}_pf_exponent = np.interp(tf.T9, xp={nuc}_temp_array, fp=np.log10({nuc}_pf_array))\n"
                    fstring += f"    {nuc}_pf = 10.0**{nuc}_pf_exponent\n"
                else:
                    fstring += f"    # setting {nuc} partition function to 1.0 by default, independent of T\n"
                    fstring += f"    {nuc}_pf = 1.0\n"
                fstring += "\n"

            fstring += "    "
            fstring += "z_r = "
            fstring += "*".join([f"{nucr}_pf" for nucr in self.rate.reactants])

            fstring += "\n"
            fstring += "    "
            fstring += "z_p = "
            fstring += "*".join([f"{nucp}_pf" for nucp in self.rate.products])

            fstring += "\n"
            fstring += f"    rate_eval.{self.fname} *= z_r/z_p\n"

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

        self._warn_about_missing_pf_tables()

        extra_args = ["[[maybe_unused]] part_fun::pf_cache_t& pf_cache", *extra_args]
        fstring = super().function_string_cxx(dtype=dtype,
                                              specifiers=specifiers,
                                              leave_open=True,
                                              extra_args=extra_args)

        # right now we have rate and drate_dT without the partition
        # function now the partition function corrections

        if self.use_pf:

            fstring += "\n"
            for nuc in set(self.rate.reactants + self.rate.products):
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
            fstring += " * ".join([f"{nucr}_pf" for nucr in self.rate.reactants])
            fstring += ";\n"

            fstring += f"    {dtype} z_p = "
            fstring += " * ".join([f"{nucp}_pf" for nucp in self.rate.products])
            fstring += ";\n\n"

            # now the derivatives, via chain rule
            chain_terms = []
            for n in self.rate.reactants:
                chain_terms.append(" * ".join([f"{nucr}_pf" for nucr in self.rate.reactants if nucr != n] + [f"d{n}_pf_dT"]))

            fstring += f"    {dtype} dz_r_dT = "
            fstring += " + ".join(chain_terms)
            fstring += ";\n"

            chain_terms = []
            for n in self.rate.products:
                chain_terms.append(" * ".join([f"{nucp}_pf" for nucp in self.rate.products if nucp != n] + [f"d{n}_pf_dT"]))

            fstring += f"    {dtype} dz_p_dT = "
            fstring += " + ".join(chain_terms)
            fstring += ";\n\n"

            fstring += f"    {dtype} dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);\n\n"

            # final terms

            fstring += "    drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);\n"
            fstring += "    rate *= z_r/z_p;\n\n"

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

        react_counts = Counter(self.rate.reactants)
        prod_counts = Counter(self.rate.products)

        reactant_factor = 1.0
        for nuc in set(self.rate.reactants):
            reactant_factor *= math.factorial(react_counts[nuc])

        product_factor = 1.0
        for nuc in set(self.rate.products):
            product_factor *= math.factorial(prod_counts[nuc])

        return (reactant_factor, product_factor)
