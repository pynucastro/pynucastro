import importlib
import os
import sys

import pytest

import pynucastro as pyna


class TestPythonPartitionNetwork:
    @pytest.fixture(scope="class")
    def pynet(self, reaclib_library):

        fwd_reactions = reaclib_library.derived_forward()

        nuclei = ["p", "he4", "fe52", "ni56", "co55"]

        fwd_rates_lib = fwd_reactions.linking_nuclei(nuclist=nuclei,
                                                     with_reverse=False)

        derived = []
        for r in fwd_rates_lib.get_rates():
            d = pyna.DerivedRate(rate=r, compute_Q=False, use_pf=True)
            derived.append(d)

        der_rates_lib = pyna.Library(rates=derived)
        full_lib = fwd_rates_lib + der_rates_lib

        return pyna.PythonNetwork(libraries=[full_lib])

    def test_partition_rates(self, pynet):
        """test the rate evaluation with partition functions from the
        python network"""
        pynet.write_network("der_net.py")
        der_net = importlib.import_module("der_net")

        T = 5.e9
        tf = pyna.Tfactors(T)

        rate_eval = der_net.RateEval()

        der_net.p_co55__he4_fe52__derived(rate_eval, tf)
        der_net.ni56__p_co55__derived(rate_eval, tf)

        assert rate_eval.p_co55__he4_fe52__derived == pytest.approx(4.570999237208017, rel=1.e-10)

        assert rate_eval.ni56__p_co55__derived == pytest.approx(23790871.179938074, rel=1.e-10)

        T = 9.e9
        tf = pyna.Tfactors(T)

        rate_eval = der_net.RateEval()

        der_net.p_co55__he4_fe52__derived(rate_eval, tf)
        der_net.ni56__p_co55__derived(rate_eval, tf)

        assert rate_eval.p_co55__he4_fe52__derived == pytest.approx(15485.753590182012, rel=1.e-10)
        assert rate_eval.ni56__p_co55__derived == pytest.approx(428973340937.6744, rel=1.e-10)

        # clean up generated files if the test passed
        os.remove("der_net.py")
        # remove imported module from cache
        del der_net
        del sys.modules["der_net"]
