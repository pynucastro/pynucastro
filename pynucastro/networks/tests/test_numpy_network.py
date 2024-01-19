import pytest
from numpy.testing import assert_allclose

import pynucastro as pyna


class TestNumpyNetwork:
    """Make sure the vectorized methods give the same results."""

    @pytest.fixture(scope="class")
    def net(self, reaclib_library):
        rate_names = ["c12(p,g)n13",
                      "c13(p,g)n14",
                      "n13(,)c13",
                      "n13(p,g)o14",
                      "n14(p,g)o15",
                      "n15(p,a)c12",
                      "o14(,)n14",
                      "o15(,)n15",
                      "he4(aa,g)c12"]
        rates = reaclib_library.get_rate_by_name(rate_names)
        net = pyna.NumpyNetwork(rates=rates)
        return net

    @pytest.fixture(scope="class")
    def comp(self, net):
        c = pyna.Composition(net.unique_nuclei)
        c.set_solar_like()
        return c

    @pytest.fixture(scope="class")
    def rho(self):
        return 1e5

    @pytest.fixture(scope="class")
    def temp(self):
        return 1e8

    def test_yfac_arr(self, net, comp):
        expected = [0.0001666666666666666, 0.0001538461538461538,
                    0.00021978021978021975, 0.0001538461538461538,
                    0.00014285714285714281, 0.00013333333333333329,
                    0.0002040816326530612, 0.00019047619047619045,
                    0.00034299999999999966]

        net.clear_arrays()
        net.update_yfac_arr(comp)

        assert_allclose(net.yfac, expected, rtol=1e-10, atol=1e-100)

    def test_prefac_arr(self, net, rho, comp):
        expected = [rho, rho, 1.0, rho, rho, rho, 1.0, 1.0, rho*rho / 6]

        net.clear_arrays()
        net.update_prefac_arr(rho, comp)

        assert_allclose(net.prefac, expected, rtol=1e-10, atol=1e-100)

    def test_evaluate_rates_arr(self, net, rho, comp, temp):
        expected = []
        rv = net.evaluate_rates(rho=rho, T=temp, composition=comp)
        expected = [rv[r] for r in net.rates]

        net.clear_arrays()
        net.update_yfac_arr(comp)
        with pytest.raises(Exception):
            net.evaluate_rates_arr(temp)
        net.update_prefac_arr(rho, comp)
        rates_arr = net.evaluate_rates_arr(temp)

        assert_allclose(rates_arr, expected, rtol=1e-10, atol=1e-100)

    def test_evaluate_ydots_arr(self, net, rho, comp, temp):
        ydots = net.evaluate_ydots(rho=rho, T=temp, composition=comp)
        expected = [ydots[nuc] for nuc in net.unique_nuclei]

        net.clear_arrays()
        with pytest.raises(Exception):
            net.evaluate_ydots_arr(temp)
        net.update_yfac_arr(comp)
        net.update_prefac_arr(rho, comp)
        ydots_arr = net.evaluate_ydots_arr(temp)

        assert_allclose(ydots_arr, expected, rtol=1e-10, atol=1e-100)

    def test_evaluate_activity_arr(self, net, rho, comp, temp):
        activity = net.evaluate_activity(rho=rho, T=temp, composition=comp)
        expected = [activity[nuc] for nuc in net.unique_nuclei]

        net.clear_arrays()
        with pytest.raises(Exception):
            net.evaluate_activity_arr(temp)
        net.update_yfac_arr(comp)
        net.update_prefac_arr(rho, comp)
        activity_arr = net.evaluate_activity_arr(temp)

        assert_allclose(activity_arr, expected, rtol=1e-10, atol=1e-100)
