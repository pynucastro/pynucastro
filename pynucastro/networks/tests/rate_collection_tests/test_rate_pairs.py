# unit tests for RateCollection.get_rate_pairs

from pynucastro.networks import RateCollection
from pynucastro.rates import Rate


class TestRatePairs:

    def test_same_mass(self):

        # we should pair these rates even though the ordering of "n"
        # and "p" switch between the two
        forward = Rate(reactants=["n", "p"], products=["d"])
        reverse = Rate(reactants=["d"], products=["p", "n"])

        net = RateCollection(rates=[forward, reverse])

        rate_pairs = net.get_rate_pairs()

        assert len(rate_pairs) == 1
