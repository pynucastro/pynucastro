from pynucastro.networks import common_networks


class TestCommon:

    def test_aprox13(self):

        net = common_networks.aprox13()

        assert len(net.unique_nuclei) == 13
        assert len(net.approx_nuclei) == 9

        assert len(net.reaclib_rates) == 37
        assert len(net.approx_rates) == 30
        assert len(net.derived_rates) == 37

    def test_mesa_basic(self):

        net = common_networks.mesa_basic()

        assert len(net.unique_nuclei) == 8
        assert len(net.approx_nuclei) == 3

        assert len(net.reaclib_rates) == 19
        assert len(net.approx_rates) == 3
        assert len(net.derived_rates) == 9
        assert len(net.branched_rates) == 2
        assert len(net.modified_rates) == 6
