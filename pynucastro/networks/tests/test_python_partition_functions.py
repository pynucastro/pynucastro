import importlib
import sys
from pathlib import Path

import pytest

import pynucastro as pyna


class TestPythonPartitionNetwork:
    @pytest.fixture(scope="class")
    def pynet(self, reaclib_library):

        fwd_reactions = reaclib_library.forward_for_detailed_balance()

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

        der_net.p_Co55__He4_Fe52__derived(rate_eval, tf)
        der_net.Ni56__p_Co55__derived(rate_eval, tf)

        assert rate_eval.p_Co55__He4_Fe52__derived == pytest.approx(4.5709992373566735, rel=1.e-10)

        assert rate_eval.Ni56__p_Co55__derived == pytest.approx(23790871.142872408, rel=1.e-10)

        T = 9.e9
        tf = pyna.Tfactors(T)

        rate_eval = der_net.RateEval()

        der_net.p_Co55__He4_Fe52__derived(rate_eval, tf)
        der_net.Ni56__p_Co55__derived(rate_eval, tf)

        assert rate_eval.p_Co55__He4_Fe52__derived == pytest.approx(15485.753590461827, rel=1.e-10)
        assert rate_eval.Ni56__p_Co55__derived == pytest.approx(428973340215.9292, rel=1.e-10)

        # clean up generated files if the test passed
        Path("der_net.py").unlink()
        # remove imported module from cache
        del der_net
        del sys.modules["der_net"]
