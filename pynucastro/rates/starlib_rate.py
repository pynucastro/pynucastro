"""Classes and methods for describing a reaction rate that is
tabulated in terms of and temperature.along with corresponding
uncertainties as per the StarLib Library

"""

import numpy as np

from pynucastro.rates.temperature_tabular_rate import TemperatureTabularRate


class StarLibRate(TemperatureTabularRate):
    """A rate whose temperature dependence and factor uncertainty are tabulated.
    Upon creation of an instance, rates at all tabulated log(T9) values are sampled.

    Parameters
    ----------
    log_t9_data : numpy.ndarray
        The temperature (in log(T / 1.e9 K)) where we tabulate the rate
    log_rate_data : numpy.ndarray
        The tabulated log(rate) data, N_A <σv>
    sigma_data : numpy.ndarray
        The tabulated log(uncertainty) data
    seed: int
        Seed to pass to rng for rate sampling
    """

    def __init__(self, log_t9_data, log_rate_data, sigma_data,
                 labelprops, seed=None, label="starlib", **kwargs):

        self.sigma_data = sigma_data
        self.log_median_rates = log_rate_data

        #Ensure same number of data points in rate and sigma
        assert (len(self.log_median_rates) == len(self.sigma_data))

        sampled_rates = self.sample_rates(seed)

        #Read in labelprops and call super
        assert isinstance(labelprops, str)
        assert len(labelprops) == 5
        self.labelprops = labelprops

        rate_source = labelprops[0:4].strip()
        #Check for electron capture
        weak_type = ''
        if rate_source == 'ec':
            weak_type = "electron_capture"

        super().__init__(log_t9_data, sampled_rates, label=label,
                         rate_source=rate_source, weak_type=weak_type,
                         **kwargs)

        #Set relevant flags
        self.weak = labelprops[4] == 'w' or rate_source == 'ec'
        self.derived_from_inverse = labelprops[4] == 'v'

    def sample_rates(self, seed=None):
        """Sample rate values as median_rate + N(0,1)*sigma for each of the
        60 entries in the data for a Starlib rate.

        Parameters
        ----------
        seed : int
            Seed for the rng necessary to sample rates. If seed is none, the
            method returns median rates.
        """

        #When no seed is provided, rates are median values
        sampled = self.log_median_rates.copy()
        #Otherwise sample using normal distribution
        if seed is not None:
            rng = np.random.default_rng(seed=seed)
            sampled = (self.log_median_rates +
                       rng.normal(size=len(self.log_median_rates))*self.sigma_data)
        return sampled

    def __eq__(self, other):
        """Determine whether two Rate objects are equal.  They are
        equal if they contain identical reactants and products.

        """

        if not isinstance(other, StarLibRate):
            return False

        return self.reactants == other.reactants and self.products == other.products
