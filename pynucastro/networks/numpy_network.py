"""Classes and methods for a network that uses NumPy arrays to cache
some intermediate data.

"""

import numpy as np

from pynucastro.networks.rate_collection import RateCollection
from pynucastro.rates import Tfactors


class NumpyNetwork(RateCollection):
    """A network that uses numpy arrays to evaluate rates more efficiently.

    Attributes
    ----------
    yfac : numpy.ndarray
        Array storing the molar fraction component of each rate (Y of each
        reactant raised to the appropriate power).  Depends on composition
        only.
    prefac : numpy.ndarray
       Array storing the prefactor for each rate, which includes both the
       statistical prefactor and mass density raised to the corresponding
       density exponent.  Depends on composition and density.

    """

    def __init__(self, rate_files=None, libraries=None, rates=None,
                 inert_nuclei=None,
                 symmetric_screening=False, do_screening=True):
        super().__init__(rate_files, libraries, rates, inert_nuclei,
                         symmetric_screening, do_screening)

        # cached values for vectorized evaluation
        self._nuc_prod_count = None
        self._nuc_cons_count = None
        self._nuc_used = None
        self._coef_arr = None
        self._coef_mask = None
        self.prefac = None
        self.yfac = None

    def _build_collection(self):
        super()._build_collection()
        # clear the cached arrays after changing any of the rates
        self.clear_arrays()

    @property
    def nuc_prod_count(self):
        """Return an array storing the count of each nucleus in rates
        producing that nucleus, with shape ``(number_of_species,
        number_of_rates)``

        """

        if self._nuc_prod_count is None:
            self._calc_count_matrices()
        return self._nuc_prod_count

    @property
    def nuc_cons_count(self):
        """Return an array storing the count of each nucleus in rates
        consuming that nucleus, with shape ``(number_of_species,
        number_of_rates)``.

        """
        if self._nuc_cons_count is None:
            self._calc_count_matrices()
        return self._nuc_cons_count

    @property
    def nuc_used(self):
        """Return a boolean matrix of whether the nucleus is involved
        in the reaction or not, with shape ``(number_of_species,
        number_of_rates)``

        """
        if self._nuc_used is None:
            self._calc_count_matrices()
        return self._nuc_used

    def _calc_count_matrices(self):
        """Compute and store 3 count matrices that can be used for
        vectorized rate calculations.

        """

        # Rate -> index mapping
        r_map = {}
        for i, r in enumerate(self.rates):
            r_map[r] = i

        N_species = len(self.unique_nuclei)
        N_rates = len(self.rates)

        # Counts for reactions producing nucleus
        self._nuc_prod_count = np.zeros((N_species, N_rates), dtype=np.int32)
        # Counts for reactions consuming nucleus
        self._nuc_cons_count = np.zeros((N_species, N_rates), dtype=np.int32)

        for i, n in enumerate(self.unique_nuclei):

            for r in self.nuclei_produced[n]:
                self._nuc_prod_count[i, r_map[r]] = r.products.count(n)

            for r in self.nuclei_consumed[n]:
                self._nuc_cons_count[i, r_map[r]] = r.reactants.count(n)

        # Whether the nucleus is involved in the reaction or not
        self._nuc_used = np.logical_or(self._nuc_prod_count, self._nuc_cons_count).T

    @property
    def coef_arr(self):
        """Reaclib rate coefficient array, with shape
        ``(number_of_rates, number_of_sets, 7)``

        """
        if self._coef_arr is None:
            self._update_rate_coef_arr()
        return self._coef_arr

    @property
    def coef_mask(self):
        """Boolean mask array determining how many sets to include in
        the final rate evaluation, with shape ``(number_of_rates,
        number_of_sets)``.

        """
        if self._coef_mask is None:
            self._update_rate_coef_arr()
        return self._coef_mask

    def _update_rate_coef_arr(self):
        """Update the :attr:`.coef_arr` and :attr:`.coef_mask` arrays."""

        # coef arr can be precomputed if evaluate_rates_arr is called multiple times
        N_sets = max(len(r.sets) for r in self.rates)

        coef_arr = np.zeros((len(self.rates), N_sets, 7), dtype=np.float64)
        coef_mask = np.zeros((len(self.rates), N_sets), dtype=np.bool_)

        for i, r in enumerate(self.rates):
            for j, s in enumerate(r.sets):
                coef_arr[i, j, :] = s.a
                coef_mask[i, j] = True

        self._coef_arr = coef_arr
        self._coef_mask = coef_mask

    def update_yfac_arr(self, composition):
        """Calculate and store molar fraction component of each rate
        (Y of each reactant raised to the appropriate power). The
        results are stored in the :attr:`.yfac` array.

        Parameters
        ----------
        composition : Composition
            the composition we are considering

        """

        # yfac must be evaluated each time composition changes, probably pretty cheap
        yfac = np.ones((len(self.rates), len(self.unique_nuclei)), dtype=np.float64)
        ys = np.array(list(composition.get_molar().values()), dtype=np.float64)

        yfac *= ys**self.nuc_cons_count.T
        yfac = np.prod(yfac, axis=1)
        self.yfac = yfac

    def update_prefac_arr(self, rho, composition):
        """Calculate and store rate prefactors, which include both
        statistical prefactors and mass density raised to the
        corresponding density exponents. The results are stored in the
        :attr:`.prefac` array.

        Parameters
        ----------
        rho : float
            mass density
        composition : Composition
            composition to evaluate rates with

        """

        y_e = composition.ye
        prefac = np.zeros(len(self.rates))
        for i, r in enumerate(self.rates):
            if r.label == "tabular":
                raise ValueError('Tabular rates are not supported in vectorized rate calculations.')
            prefac[i] = r.prefactor * rho**r.dens_exp
            if r.weak_type == 'electron_capture':
                prefac[i] *= y_e

        self.prefac = prefac

    def evaluate_rates_arr(self, T):
        """Evaluate the rates in the network for a specific
        temperature, assuming necessary precalculations have been
        carried out (calling the methods :meth:`.update_yfac_arr` and
        :meth:`.update_prefac_arr` to set the composition and
        density).

        This performs a vectorized calculation, and returns an array ordered by
        the rates in the ``rates`` member variable. This does not support
        tabular rates.

        See :meth:`.evaluate_rates` for the non-vectorized version. Relative
        performance between the two varies based on the setup. See
        :meth:`clear_arrays` for freeing memory post calculation.

        Parameters
        ----------
        T : float
            temperature

        Returns
        -------
        numpy.ndarray

        """

        # T9 arr only needs to be evaluated when T changes
        T9_arr = Tfactors(T).array[None, None, :]

        rvals = self.prefac*self.yfac*np.sum(np.exp(np.sum(self.coef_arr*T9_arr, axis=2))*self.coef_mask, axis=1)

        return rvals

    def evaluate_ydots_arr(self, T):
        """Evaluate net rate of change of molar abundance for each nucleus in the
        network for a specific temperature, assuming necessary precalculations
        have been carried out (calling the methods :meth:`.update_yfac_arr` and
        :meth:`.update_prefac_arr` to set the composition and density).

        This performs a vectorized calculation, and returns an array
        ordered by the nuclei in the ``unique_nuclei`` member
        variable. This does not support tabular rates.

        See :meth:`.evaluate_ydots` for the non-vectorized version. Relative
        performance between the two varies based on the setup. See
        :meth:`.clear_arrays` for freeing memory post calculation.

        Parameters
        ----------
        T : float
            temperature

        Returns
        -------
        numpy.ndarray

        """

        rvals_arr = self.evaluate_rates_arr(T)

        p_A = np.sum(self.nuc_prod_count*rvals_arr, axis=1)
        c_A = np.sum(self.nuc_cons_count*rvals_arr, axis=1)

        return p_A - c_A

    def evaluate_activity_arr(self, T):
        """Sum over all of the terms contributing to dY/dt for a
        specific temperature, neglecting sign, assuming necessary
        precalculations have been carried out (calling the methods
        :meth:`.update_yfac_arr` and :meth:`.update_prefac_arr` to set
        the composition and density).

        This performs a vectorized calculation, and returns an array ordered by
        the nuclei in the ``unique_nuclei`` member variable. This does not
        support tabular rates.

        See :meth:`.evaluate_activity` for the non-vectorized version. Relative
        performance between the two varies based on the setup. See
        :meth:`.clear_arrays` for freeing memory post calculation.

        Parameters
        ----------
        T : float
            temperature

        Returns
        -------
        numpy.ndarray

        """

        rvals_arr = self.evaluate_rates_arr(T)

        p_A = np.sum(self.nuc_prod_count*rvals_arr, axis=1)
        c_A = np.sum(self.nuc_cons_count*rvals_arr, axis=1)

        return p_A + c_A

    def clear_arrays(self):
        """Clear all cached arrays stored by the
        :meth:`.update_yfac_arr` and :meth:`.update_prefac_arr` member
        functions, freeing up memory.

        """
        self._nuc_prod_count = None
        self._nuc_cons_count = None
        self._nuc_used = None
        self._coef_arr = None
        self._coef_mask = None
        self.prefac = None
        self.yfac = None
