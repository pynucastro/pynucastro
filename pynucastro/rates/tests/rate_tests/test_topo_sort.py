# Test the topological sorting of rates by setting up a network with
# different rates types that carry dependencies.

import copy
import random

import pytest

from pynucastro.networks import RateCollection
from pynucastro.nucdata import Nucleus
from pynucastro.rates import BranchedRate, ModifiedRate, aprox_family_rates
from pynucastro.sort_utils import CircularRateDependency, topo_sort


class TestTopoSort:
    @pytest.fixture(scope="class")
    @classmethod
    def net(cls, reaclib_library):

        # do the simple CNO approximation

        rc12pg = reaclib_library.get_rate_by_name("c12(p,g)n13")
        rn14pg = reaclib_library.get_rate_by_name("n14(p,g)o15")
        rn15pa = reaclib_library.get_rate_by_name("n15(p,a)c12")
        rn15pg = reaclib_library.get_rate_by_name("n15(p,g)o16")
        ro16pg = reaclib_library.get_rate_by_name("o16(p,g)f17")

        rc12_2p_n14 = ModifiedRate(rc12pg,
                                   new_products=[Nucleus("n14")],
                                   stoichiometry={Nucleus("p"): 2})

        rn14_2p_c12 = BranchedRate(rn14pg,
                                   primary_branch=rn15pa,
                                   other_branch=rn15pg,
                                   stoichiometry={Nucleus("p"): 2},
                                   description="N14(p,g)O15(,e+nu)N15(p,a)C12")

        rn14_2p_o16 = BranchedRate(rn14pg,
                                   primary_branch=rn15pg,
                                   other_branch=rn15pa,
                                   stoichiometry={Nucleus("p"): 2},
                                   description="N14(p,g)O15(,e+nu)N15(p,g)O16")

        ro16_2p_n14_a = ModifiedRate(ro16pg,
                                     new_products=[Nucleus("n14"),
                                                   Nucleus("he4")],
                                     stoichiometry={Nucleus("p"): 2})

        # add add in the C12(a,g)O16 approximation

        rc12ag, _ = aprox_family_rates.make_ap_pg_rates(reaclib_library, "c12", "o16")

        net = RateCollection(rates=[rc12_2p_n14,
                                    rn14_2p_c12,
                                    rn14_2p_o16,
                                    ro16_2p_n14_a,
                                    rc12ag])

        return net

    def test_dependencies(self, net):

        all_rates = list(net.all_rates)

        # randomly sort the rates in the network and then do a
        # topological sort and check to make sure all dependencies
        # come before the rates that depend on them.

        random.seed(1234)
        random.shuffle(all_rates)

        sorted_rates = topo_sort(all_rates)

        # now check that there for each rate in the list
        # their dependencies come before them

        for i, r in enumerate(sorted_rates):
            if crates := r.get_child_rates():
                for cr in crates:
                    assert cr in sorted_rates[:i]

    def test_circular_dependency(self, reaclib_library, net):

        all_rates = list(net.all_rates)

        # we will fake a circular dependency by hacking get_child_rates in an instance
        rn14pg = reaclib_library.get_rate_by_name("n14(p,g)o15")
        all_rates.remove(rn14pg)

        rn14pg_new = copy.copy(rn14pg)

        rn14_2p_c12 = [r for r in all_rates
                       if r.reactants == [Nucleus("p"), Nucleus("n14")] and
                          r.products == [Nucleus("he4"), Nucleus("c12")]][0]

        rn14pg_new.get_child_rates = lambda: [rn14_2p_c12]

        all_rates.append(rn14pg_new)

        # randomly sort the rates in the network and then do a
        # topological sort and check to make sure all dependencies
        # come before the rates that depend on them.

        random.seed(1234)
        random.shuffle(all_rates)

        # this will fail with a circular dependency

        with pytest.raises(CircularRateDependency):
            sorted_rates = topo_sort(all_rates)
