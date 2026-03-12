"""Classes and methods for describing a reaction rate that is
tabulated in terms of and temperature.along with corresponding
uncertainties as per the StarLib Library

"""

import numpy as np

from pynucastro.rates.temperature_tabular_rate import TemperatureTabularRate, TempTableInterpolator


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
                 labelprops, rng=None, label="starlib", **kwargs):

        #Read in labelprops and call super
        assert isinstance(labelprops, str)
        assert len(labelprops) == 5
        self.labelprops = labelprops

        rate_source = labelprops[0:4].strip()
        #Check for electron capture
        weak_type = ''
        if rate_source == 'ec':
            weak_type = "electron_capture"

        super().__init__(log_t9_data, log_rate_data, label=label,
                         rate_source=rate_source, weak_type=weak_type,
                         **kwargs)

        #Set relevant flags
        self.weak = labelprops[4] == 'w' or rate_source == 'ec'
        self.derived_from_inverse = labelprops[4] == 'v'

        #Store data necessary for sampling
        self.sigma_data = sigma_data
        self.log_median_rates = log_rate_data.copy()

        #Redefine the Interpolator to reflect sampled rates
        self.sample_rates(rng=rng)

        #Ensure same number of data points in rate and sigma
        assert (len(self.log_median_rates) == len(self.sigma_data))

    def sample_rates(self, rng=None):
        """Sample rate values as median_rate + N(0,1)*sigma for each of the
        60 entries in the data for a Starlib rate, given a non-empty seed. Use
        sampled rates to overwrite interpolator

        Parameters
        ----------
        rng : np.random.default_rng
            An rng that draws the gaussian deviate required for rate sampling
        """

        #When no rng is provided, rates are median values
        sampled = self.log_median_rates.copy()

        #Otherwise sample using normal distribution
        if rng is not None:
            p = rng.normal()
            sampled += p*self.sigma_data

        #Rewrite interpolator
        self.log_rate_data = sampled
        self.interpolator = TempTableInterpolator(self.log_t9_data, sampled)

    def __eq__(self, other):
        """Determine whether two Rate objects are equal.  They are
        equal if they contain identical reactants and products.

        """

        if not isinstance(other, StarLibRate):
            return False

        return self.reactants == other.reactants and self.products == other.products
