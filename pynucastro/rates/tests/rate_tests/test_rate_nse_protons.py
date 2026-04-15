# test the rate swap_protons method

import copy

import pynucastro as pyna


class TestRate:

    def test_swap_protons(self, reaclib_library):

        _mn51pg = reaclib_library.get_rate_by_name("mn51(p,g)fe52")
        mn51pg = copy.deepcopy(_mn51pg)

        p = pyna.Nucleus("p")
        p_nse = pyna.Nucleus("p_nse")

        assert p in mn51pg.reactants
        assert p_nse not in mn51pg.reactants

        mn51pg.swap_protons()

        assert p not in mn51pg.reactants
        assert p_nse in mn51pg.reactants
