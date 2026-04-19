"""Classes and methods for describing a reaction rate that is
tabulated in terms of and temperature.along with corresponding
uncertainties as per the StarLib Library

"""

from pynucastro.rates.temperature_tabular_rate import (TemperatureTabularRate,
                                                       TempTableInterpolator)


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
    rng : numpy.random.Generator
        An rng for rate sampling

    Attributes
    ----------
    log_median_rate : numpy.ndarray
        A backup of the median rate, stored so we can resample

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
        rng : numpy.random.Generator
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

        return (self.reactants == other.reactants and
                self.products == other.products and
                self.weak_type == other.weak_type)

    def __hash__(self):
        return hash(self.__repr__())

    def function_string_cxx(self, dtype="double", specifiers="inline",
                            leave_open=False, extra_args=None):
        """Return a string containing the C++ function that computes
        the rate, taking into account the uncertainty.

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

        # pylint: disable=duplicate-code
        if extra_args is None:
            extra_args = ()

        args = ["const tf_t& tfactors",
                f"const {dtype} log_scor", f"const {dtype} dlog_scor_dT",
                f"{dtype}& rate", f"{dtype}& drate_dT", *extra_args]
        fstring = ""
        fstring += "template <int do_T_derivatives>\n"
        fstring += f"{specifiers}\n"
        fstring += f"void rate_{self.fname}({', '.join(args)}) {{\n\n"
        fstring += f"    // {self.rid}\n\n"
        # pylint: enable=duplicate-code
        fstring += "    // our rate is exp(μ + pσ + h)\n"
        fstring += "    // where μ = median rate, p = Gaussian random #,\n"
        fstring += "    //       σ = uncertainty, h = screening potential\n"
        fstring += "    auto [_mu, _dmu_dlogT9] = interp_net::cubic_interp_uneven<do_T_derivatives>(\n"
        fstring += "                                          tfactors.lnT9,\n"
        fstring += f"                                          {self.fname}_data::log_t9,\n"
        fstring += f"                                          {self.fname}_data::log_rate);\n"
        fstring += f"   auto p = Rates::get_p_random<k_{self.fname}>();\n"
        fstring += "    auto [_sigma, _dsigma_dlogT9] = interp_net::cubic_interp_uneven<do_T_derivatives>(\n"
        fstring += "                                                 tfactors.lnT9,\n"
        fstring += f"                                                 {self.fname}_data::log_t9,\n"
        fstring += f"                                                 {self.fname}_data::sigma_rate);\n"
        fstring += "    rate = std::exp(_mu + p * _sigma + log_scor);\n"
        fstring += "    // we found dlog(rate)/dlog(T9)\n"
        fstring += "    if constexpr (do_T_derivatives) {\n"
        fstring += f"        {dtype} dlog_rate_dT = tfactors.T9i * 1.e-9_rt * (_dmu_dlogT9 + p * _dsigma_dlogT9) + dlog_scor_dT\n;"
        fstring += "        drate_dT = rate * dlog_rate_dT;\n"
        fstring += "    }\n"

        if not leave_open:
            fstring += "}\n\n"

        return fstring
