# this will create a network using a ModifiedRate.

import copy

import numpy as np
import pytest
from pytest import approx

import pynucastro as pyna


class TestModifiedRate:
    @pytest.fixture(scope="class")
    @classmethod
    def new_net(cls, reaclib_library):
        # create a network that uses a ModifiedRate

        lib = reaclib_library.linking_nuclei(["he4", "c12", "o16",
                                              "ne20", "mg24"],
                                             with_reverse=False)

        _c12c12_other = reaclib_library.get_rate_by_name("c12(c12,n)mg23")
        c12c12_other = copy.deepcopy(_c12c12_other)
        c12c12_new = pyna.ModifiedRate(c12c12_other,
                                       new_products=["mg24"])
        lib.add_rate(c12c12_new)

        return pyna.PythonNetwork(libraries=[lib])

    def test_hidden_rates(self, new_net):

        assert len(new_net.get_hidden_rates()) == 1

    def test_all_rates(self, new_net):

        assert len(new_net.all_rates) == 8

    def test_module(self, new_net):
        """write the new network to a file and import it and then
        compare the ydots that the rhs() function gives to the
        evaluate_ydots from the PythonNetwork object"""

        comp = pyna.Composition(new_net.unique_nuclei)
        comp.set_equal()

        rho = 1.e8
        T = 2.e9

        new_net.write_network("test_modified_net.py")
        import test_modified_net as mn  # pylint: disable=import-outside-toplevel,import-error  # noqa: PLC0415

        Y = np.asarray(list(comp.get_molar().values()))

        module_ydots = mn.rhs(0.0, Y, rho, T)

        state = pyna.ThermoState(rho=rho, T=T, comp=comp)
        net_ydots = new_net.evaluate_ydots(state)

        for n, k in enumerate(net_ydots):
            assert net_ydots[k] == approx(module_ydots[n], rel=1.e-11, abs=1.e-14)

    def test_Q_value(self, reaclib_library):
        r = reaclib_library.get_rate_by_name("c12(a,g)o16")
        mr_1 = pyna.ModifiedRate(r, new_reactants=["he4", "c12"],
                                 new_products=["ne20"], stoichiometry={pyna.Nucleus("he4"): 2})
        mr_2 = pyna.ModifiedRate(r, new_reactants=["he4", "he4", "c12"], new_products=["ne20"])

        assert mr_1.Q == mr_2.Q
