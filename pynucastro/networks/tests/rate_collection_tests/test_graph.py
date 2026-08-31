import warnings

import numpy as np
import pytest
from pytest import approx

from pynucastro.networks import network_helper
from pynucastro.nucdata import Composition, Nucleus
from pynucastro.rates import ThermoState

# tests of the networkx graph function




class TestNetworkXGraph:

    @pytest.fixture(scope="class")
    @classmethod
    def net(cls):
        nuclei = ["p", "he4",
                  "c12", "c13",
                  "n13", "n14", "n15",
                  "o14", "o15", "o16", "o17", "o18",
                  "f17", "f18", "f19",
                  "ne18", "ne19", "ne20"]

        return network_helper(nuclei, tabular_ordering=["ffn", "oda"])

    @pytest.fixture(scope="class")
    @classmethod
    def state(cls, net):
        rho = 1.e4
        T = 5.e8
        comp = Composition(net.unique_nuclei, init="solar")
        return ThermoState(rho=rho, T=T, comp=comp)

    def test_probability(self, net, state):


        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=UserWarning)
            rates = net.evaluate_rates(state)

        G = net.create_network_graph(net.unique_nuclei,
                                     rate_ydots=rates,
                                     use_branching_ratios=True)

        # check that the branching probabilities leaving a nucleus sum
        # to 1.

        for n in net.unique_nuclei:
            # we hide p and alpha links by default, so the
            # probabilities won't add to 1
            if n in [Nucleus("p"), Nucleus("he4")]:
                continue

            ptot = 0.0
            for u, _, data in G.edges(data=True):
                if u == n:
                    ptot += np.exp(-data["weight"])

            assert ptot == approx(1.0)

