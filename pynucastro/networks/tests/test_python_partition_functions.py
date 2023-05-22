import pytest

import pynucastro as pyna


class TestPythonPartitionNetwork:
    @pytest.fixture(scope="class")
    def fn(self, reaclib_library):

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

        pynet = pyna.PythonNetwork(libraries=[full_lib])
        pynet.write_network("der_net.py")

    def test_partition_rates(self, fn):
        """test the rate evaluation with partition functions from the
        python network"""

        import der_net

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
