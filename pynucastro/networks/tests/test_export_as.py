import pytest

import pynucastro as pyna


class TestNetworkExport:
    @pytest.fixture(scope="class")
    def net(self):
        _net =  pyna.network_helper(["p", "he4", "c12", "o16",
                                     "ne20", "na23", "mg24"])

        return pyna.PythonNetwork(rates=_net.get_rates(),
                                  inert_nuclei=["ni56"])

    def test_convert_cxx(self, net):
        assert isinstance(net, pyna.PythonNetwork)

        new_net = net.export_as(pyna.SimpleCxxNetwork)

        assert isinstance(new_net, pyna.SimpleCxxNetwork)

        assert len(net.unique_nuclei) == len(new_net.unique_nuclei)
        assert len(net.inert_nuclei) == len(new_net.inert_nuclei)

        assert len(net.reaclib_rates) == len(new_net.reaclib_rates)
        assert len(net.derived_rates) == len(new_net.derived_rates)
