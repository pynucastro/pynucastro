# unit tests for rates

from pynucastro.nucdata import Nucleus


class TestRate:

    def test_swap_protons(self, reaclib_library):

        mn51pg = reaclib_library.get_rate_by_name("mn51(p,g)fe52")

        p = Nucleus("p")
        p_nse = Nucleus("p_nse")

        assert p in mn51pg.reactants
        assert p_nse not in mn51pg.reactants

        mn51pg.swap_protons()

        assert p not in mn51pg.reactants
        assert p_nse in mn51pg.reactants
