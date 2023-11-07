import numpy as np

from pynucastro.networks.rate_collection import RateCollection
from pynucastro.rates import Tfactors


class NumpyNetwork(RateCollection):
    """A network that uses numpy arrays to evaluate rates more efficiently."""

    def __init__(self, rate_files=None, libraries=None, rates=None,
                 inert_nuclei=None,
                 symmetric_screening=False, do_screening=True):
        super().__init__(rate_files, libraries, rates, inert_nuclei,
                         symmetric_screening, do_screening)

        # cached values for vectorized evaluation
        self.nuc_prod_count = None
        self.nuc_cons_count = None
        self.nuc_used = None
        self.coef_arr = None
        self.coef_mask = None
        self.prefac = None
        self.yfac = None

    def calc_count_matrices(self):
        """
        Compute and store 3 count matrices that can be used for vectorized rate calculations.
        Each matrix has shape *number_of_species × number_of_rates*. The first matrix (*nuc_prod_count*)
        stores the count of each nucleus in rates producing that nucleus. The second
        (*nuc_cons_count*) stores the count of each nucleus in rates consuming that nucleus. The
        third (*nuc_used*) stores a Boolean matrix of whether the nucleus is involved in the reaction
        or not.
        """

        # Rate -> index mapping
        r_map = {}
        for i, r in enumerate(self.rates):
            r_map[r] = i

        N_species = len(self.unique_nuclei)
        N_rates = len(self.rates)

        # Counts for reactions producing nucleus
        self.nuc_prod_count = np.zeros((N_species, N_rates), dtype=np.int32)
        # Counts for reactions consuming nucleus
        self.nuc_cons_count = np.zeros((N_species, N_rates), dtype=np.int32)

        for i, n in enumerate(self.unique_nuclei):

            for r in self.nuclei_produced[n]:
                self.nuc_prod_count[i, r_map[r]] = r.products.count(n)

            for r in self.nuclei_consumed[n]:
                self.nuc_cons_count[i, r_map[r]] = r.reactants.count(n)

        # Whether the nucleus is involved in the reaction or not
        self.nuc_used = np.logical_or(self.nuc_prod_count, self.nuc_cons_count).T

    def update_rate_coef_arr(self):
        """
        Store Reaclib rate coefficient array, as well as a Boolean mask array determining
        how many sets to include in the final rate evaluation. The first array (*coef_arr*)
        has shape *number_of_rates × number_of_species × 7*, while the second has shape
        *number_of_rates × number_of_species*.
        """

        # coef arr can be precomputed if evaluate_rates_arr is called multiple times
        N_sets = max(len(r.sets) for r in self.rates)

        coef_arr = np.zeros((len(self.rates), N_sets, 7), dtype=np.float64)
        coef_mask = np.zeros((len(self.rates), N_sets), dtype=np.bool_)

        for i, r in enumerate(self.rates):
            for j, s in enumerate(r.sets):
                coef_arr[i, j, :] = s.a
                coef_mask[i, j] = True

        self.coef_arr = coef_arr
        self.coef_mask = coef_mask

    def update_yfac_arr(self, composition):
        """
        Calculate and store molar fraction component of each rate (Y of each reactant raised to
        the appropriate power). The results are stored in an array called *yfac*. The method
        *calc_count_matrices* needs to have been called with this net beforehand.
        """

        # yfac must be evaluated each time composition changes, probably pretty cheap
        yfac = np.ones((len(self.rates), len(self.unique_nuclei)), dtype=np.float64)
        ys = np.array(list(composition.get_molar().values()), dtype=np.float64)

        yfac *= ys**self.nuc_cons_count.T
        yfac = np.prod(yfac, axis=1)
        self.yfac = yfac

    def update_prefac_arr(self, rho, composition):
        """
        Calculate and store rate prefactors, which include both statistical prefactors and
        mass density raised to the corresponding density exponents. The results are stored in
        an array called *prefac*.
        """

        y_e = composition.eval_ye()
        prefac = np.zeros(len(self.rates))
        for i, r in enumerate(self.rates):
            if r.label == "tabular":
                raise ValueError('Tabular rates are not supported in vectorized rate calculations.')
            prefac[i] = r.prefactor * rho**r.dens_exp
            if r.weak_type == 'electron_capture':
                prefac[i] *= y_e

        self.prefac = prefac

    def evaluate_rates_arr(self, T):
        """
        Evaluate the rates in the network for a specific temperature, assuming necessary precalculations
        have been carried out (calling the methods *calc_count_matrices*, *update_rate_coef_arr*,
        *update_yfac_arr*, and *update_prefac_arr*). The latter two set the composition and density.

        This performs a vectorized calculation, and returns an array ordered by the rates in the *rates*
        member variable. This does not support tabular rates.

        See *evaluate_rates* method for non-vectorized version. Relative performance between the
        two varies based on the setup. See *clear_arrays* for freeing memory post calculation.
        """

        # T9 arr only needs to be evaluated when T changes
        T9_arr = Tfactors(T).array[None, None, :]

        rvals = self.prefac*self.yfac*np.sum(np.exp(np.sum(self.coef_arr*T9_arr, axis=2))*self.coef_mask, axis=1)

        return rvals

    def evaluate_ydots_arr(self, T):
        """
        Evaluate net rate of change of molar abundance for each nucleus in the network for a
        specific temperature, assuming necessary precalculations have been carried out (calling the
        methods *calc_count_matrices*, *update_rate_coef_arr*, *update_yfac_arr*, and
        *update_prefac_arr*). The latter two set the composition and density.

        This performs a vectorized calculation, and returns an array ordered by the nuclei in the
        *unique_nuclei* member variable. This does not support tabular rates.

        See *evaluate_ydots* method for non-vectorized version. Relative performance between the
        two varies based on the setup. See *clear_arrays* for freeing memory post calculation.
        """

        rvals_arr = self.evaluate_rates_arr(T)

        p_A = np.sum(self.nuc_prod_count*rvals_arr, axis=1)
        c_A = np.sum(self.nuc_cons_count*rvals_arr, axis=1)

        return p_A - c_A

    def evaluate_activity_arr(self, T):
        """
        sum over all of the terms contributing to dY/dt for a specific temperature, neglecting sign,
        assuming necessary precalculations have been carried out (calling the methods
        *calc_count_matrices*, *update_rate_coef_arr*, *update_yfac_arr*, and *update_prefac_arr*).
        The latter two set the composition and density.

        This performs a vectorized calculation, and returns an array ordered by the nuclei in the
        *unique_nuclei* member variable. This does not support tabular rates.

        See *evaluate_activity* method for non-vectorized version. Relative performance between the
        two varies based on the setup. See *clear_arrays* for freeing memory post calculation.
        """

        rvals_arr = self.evaluate_rates_arr(T)

        p_A = np.sum(self.nuc_prod_count*rvals_arr, axis=1)
        c_A = np.sum(self.nuc_cons_count*rvals_arr, axis=1)

        return p_A + c_A

    def clear_arrays(self):
        """
        Clear all temporary variables created/set by the *calc_count_matrices*, *update_rate_coef_arr*,
        *update_yfac_arr*, and *update_prefac_arr* member functions, freeing up memory.
        """

        try:
            del self.nuc_prod_count
        except AttributeError:
            pass

        try:
            del self.nuc_cons_count
        except AttributeError:
            pass

        try:
            del self.nuc_used
        except AttributeError:
            pass

        try:
            del self.coef_arr
        except AttributeError:
            pass

        try:
            del self.coef_mask
        except AttributeError:
            pass

        try:
            del self.prefac
        except AttributeError:
            pass

        try:
            del self.yfac
        except AttributeError:
            pass
