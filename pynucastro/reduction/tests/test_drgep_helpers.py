"""Test the various helper functions in the reduction module."""

import pytest
from numpy.testing import assert_allclose

import pynucastro as pyna
import pynucastro.reduction.drgep_module as drgep


class TestDrgepHelpers:
    @pytest.fixture(scope="class")
    def net(self, reaclib_library):
        rate_names = ["c12(p,g)n13",
                      "c13(p,g)n14",
                      "n13(,)c13",
                      "n13(p,g)o14",
                      "n14(p,g)o15",
                      "n15(p,a)c12",
                      "o14(,)n14",
                      "o15(,)n15"]
        rates = reaclib_library.get_rate_by_name(rate_names)
        return pyna.NumpyNetwork(rates=rates)

    @pytest.fixture(scope="class")
    def comp(self, net):
        c = pyna.Composition(net.unique_nuclei)
        c.set_solar_like()
        return c

    def test_interaction_matrix(self, net, comp):
        expected = [
            [1.0, 0.9979577882382153, 0.9984575204946522,
             0.0014759653755533888, 0.0005550559979812046, 0.0014871557638034168,
             0.9979577882382153, 5.532374154433225e-05, 1.1190388250027898e-05],
            [1.0, 1.0, 1.0,
             0.0, 0.0, 0.0,
             1.0, 0.0, 0.0],
            [0.9994992450960084, 1.0, 0.9994992450960084,
             0.0, 0.0005007549039915752, 0.0,
             1.0, 0.0, 0.0],
            [1.0, 0.0, 0.0,
             0.9980318855335509, 0.0019681144664490413, 1.0,
             0.0, 0.0, 0.0],
            [0.8892932348638157, 0.0, 1.0,
             0.005812850321723165, 0.8834803845420925, 0.0,
             0.0, 0.1107067651361843, 0.0],
            [0.9772901230569564, 0.0, 0.0,
             0.9847562909016563, 0.0, 0.9925338321553002,
             0.0, 0.01524370909834373, 0.007466167844699782],
            [1.0, 1.0, 1.0,
             0.0, 0.0, 0.0,
             0.9999876346069883, 0.0, 1.2365393011649987e-05],
            [1.0, 0.0, 0.0,
             0.0, 1.0, 0.4129776173198909,
             0.0, 0.5870223826801091, 0.0],
            [0.906828286689306, 0.0, 0.0,
             0.0, 0.0, 0.906828286689306,
             1.0, 0.0, 0.09317171331069393],
        ]

        rvals = net.evaluate_rates(rho=1e4, T=1e8, composition=comp)
        r_AB = drgep.calc_interaction_matrix(net, rvals)
        assert_allclose(r_AB, expected, rtol=1e-10, atol=1e-100)

    def test_interaction_matrix_numpy(self, net, comp):
        rvals = net.evaluate_rates(rho=1e5, T=1e8, composition=comp)
        rvals_arr = [rvals[r] for r in net.rates]
        expected = drgep.calc_interaction_matrix(net, rvals)

        net.clear_arrays()
        r_AB = drgep.calc_interaction_matrix_numpy(net, rvals_arr)

        assert_allclose(r_AB, expected, rtol=1e-10, atol=1e-100)

    def test_adj_nuc(self, net):
        p = pyna.Nucleus("p")
        a = pyna.Nucleus("a")
        c12 = pyna.Nucleus("c12")
        c13 = pyna.Nucleus("c13")
        n13 = pyna.Nucleus("n13")
        n14 = pyna.Nucleus("n14")
        n15 = pyna.Nucleus("n15")
        o14 = pyna.Nucleus("o14")
        o15 = pyna.Nucleus("o15")
        expected = {
            p: {p, c12, n13, c13, n14, o14, o15, n15, a},
            a: {a, n15, p, c12},
            c12: {c12, p, n13, n15, a},
            c13: {c13, p, n13, n14},
            n13: {n13, c12, p, c13, o14},
            n14: {n14, c13, p, o15, o14},
            n15: {n15, p, a, c12, o15},
            o14: {o14, n13, p, n14},
            o15: {o15, n14, p, n15},
        }

        adj_nuc = drgep.get_adj_nuc(net)
        assert adj_nuc == expected
