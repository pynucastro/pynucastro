# unit tests for rates

import pynucastro as pyna


class TestRate:

    def test_swap_protons(self):

        # note: we can't use the reaclib_library fixture we've
        # setup in pytest because we will be modifying rates,
        # which will be reflected in the fixture and affect
        # later tests
        rl = pyna.ReacLibLibrary()

        mn51pg = rl.get_rate_by_name("mn51(p,g)fe52")

        p = pyna.Nucleus("p")
        p_nse = pyna.Nucleus("p_nse")

        assert p in mn51pg.reactants
        assert p_nse not in mn51pg.reactants

        mn51pg.swap_protons()

        assert p not in mn51pg.reactants
        assert p_nse in mn51pg.reactants
