import pytest

import pynucastro as pyna


class TestFullPythonNetwork:
    @pytest.fixture(scope="class")
    def fn(self, reaclib_library, tabular_library):
        rate_names = ["c12(c12,a)ne20",
                      "c12(c12,n)mg23",
                      "c12(c12,p)na23",
                      "c12(a,g)o16",
                      "n(,)p",
                      "he4(aa,g)c12"]
        rates = reaclib_library.get_rate_by_name(rate_names)

        tabular_rate_names = ["na23(,)ne23",
                              "ne23(,)na23"]
        tabular_rates = tabular_library.get_rate_by_name(tabular_rate_names)

        fn = pyna.PythonNetwork(rates=rates+tabular_rates)
        return fn

    def test_ydot_eval(self, fn):

        # these conditions were picked to make the Na23 rate from
        # ReacLib and tabulated rates comparable
        rho = 1.e8
        T = 1.2e9

        comp = pyna.Composition(fn.unique_nuclei)
        comp.set_equal()

        full_ydots = fn.evaluate_ydots(rho, T, comp)

        rl_ydots = fn.evaluate_ydots(rho, T, comp,
                                     rate_filter=lambda r: isinstance(r, pyna.rates.ReacLibRate))

        tl_ydots = fn.evaluate_ydots(rho, T, comp,
                                     rate_filter=lambda r: isinstance(r, pyna.rates.TabularRate))

        print(full_ydots)

        # these have no tabular contribution
        assert rl_ydots[pyna.Nucleus("n")] == full_ydots[pyna.Nucleus("n")]
        assert rl_ydots[pyna.Nucleus("p")] == full_ydots[pyna.Nucleus("p")]
        assert rl_ydots[pyna.Nucleus("he4")] == full_ydots[pyna.Nucleus("he4")]
        assert rl_ydots[pyna.Nucleus("c12")] == full_ydots[pyna.Nucleus("c12")]
        assert rl_ydots[pyna.Nucleus("o16")] == full_ydots[pyna.Nucleus("o16")]
        assert rl_ydots[pyna.Nucleus("mg23")] == full_ydots[pyna.Nucleus("mg23")]

        # this has no ReacLib contribution
        assert tl_ydots[pyna.Nucleus("ne23")] == full_ydots[pyna.Nucleus("ne23")]

        # this has both
        assert rl_ydots[pyna.Nucleus("na23")] + tl_ydots[pyna.Nucleus("na23")] == full_ydots[pyna.Nucleus("na23")]
