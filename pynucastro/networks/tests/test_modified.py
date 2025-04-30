# this will create 2 versions of a network with a modified rate.  the
# first will simply use Rate.modifiy_products to change the endpoints.
# The second will create a proper ModifiedRate.  We test to make sure
# that these give the same values for the rates.

import pytest
from pytest import approx

import pynucastro as pyna

import numpy as np


class TestModifiedRate:
    @pytest.fixture(scope="class")
    def original_net(self):
        # create a network and use the rate .modify_products()
        # to change the endpoints

        # note: because we are modifying a rate from ReacLibLibrary
        # directly, we cannot use the reaclib_library fixture,
        # since that change will propagate to the other test

        rl = pyna.ReacLibLibrary()

        lib = rl.linking_nuclei(["he4", "c12", "o16",
                                 "ne20", "mg24"],
                                with_reverse=False)

        c12c12_other = rl.get_rate_by_name("c12(c12,n)mg23")
        c12c12_other.modify_products(["mg24"])
        lib.add_rate(c12c12_other)

        return pyna.PythonNetwork(libraries=[lib])

    @pytest.fixture(scope="class")
    def new_net(self, reaclib_library):
        # create a network that uses a ModifiedRate

        lib = reaclib_library.linking_nuclei(["he4", "c12", "o16",
                                              "ne20", "mg24"],
                                             with_reverse=False)

        c12c12_other = reaclib_library.get_rate_by_name("c12(c12,n)mg23")
        c12c12_new = pyna.ModifiedRate(c12c12_other,
                                       new_products=["mg24"])
        lib.add_rate(c12c12_new)

        return pyna.PythonNetwork(libraries=[lib])

    def test_hidden_rates(self, original_net, new_net):

        assert len(original_net.get_hidden_rates()) == 0
        assert len(new_net.get_hidden_rates()) == 1

    def test_all_rates(self, original_net, new_net):

        assert len(original_net.all_rates) == 7
        assert len(new_net.all_rates) == 8

    def test_evaluate_rates(self, original_net, new_net):
        """test that the old and new ways of doing the modification
        give the same rates

        """

        comp = pyna.Composition(original_net.unique_nuclei)
        comp.set_equal()

        rho = 1.e8
        T = 2.e9

        onet_eval = original_net.evaluate_rates(rho=rho, T=T, composition=comp)
        nnet_eval = new_net.evaluate_rates(rho=rho, T=T, composition=comp)

        for key in onet_eval:
            assert onet_eval[key] == nnet_eval[key]

    def test_ydots(self, original_net, new_net):
        """test that the old and new ways of doing the modification
        give the same ydots

        """

        comp = pyna.Composition(original_net.unique_nuclei)
        comp.set_equal()

        rho = 1.e8
        T = 2.e9

        onet_ydot = original_net.evaluate_ydots(rho=rho, T=T, composition=comp)
        nnet_ydot = new_net.evaluate_ydots(rho=rho, T=T, composition=comp)

        for key in onet_ydot:
            assert onet_ydot[key] == nnet_ydot[key]

    def test_module(self, new_net):
        """write the new network to a file and import it and then
        compare the ydots that the rhs() function gives to the
        evaluate_ydots from the PythonNetwork object"""

        comp = pyna.Composition(new_net.unique_nuclei)
        comp.set_equal()

        rho = 1.e8
        T = 2.e9

        new_net.write_network("test_modified_net.py")
        import test_modified_net as mn

        Y = np.asarray([v for v in comp.get_molar().values()])

        module_ydots = mn.rhs(0.0, Y, rho, T)

        net_ydots = new_net.evaluate_ydots(rho=rho, T=T, composition=comp)

        for n, k in enumerate(net_ydots):
            assert net_ydots[k] == approx(module_ydots[n], rel=1.e-11, abs=1.e-14)
