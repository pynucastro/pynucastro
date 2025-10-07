import pytest

import pynucastro as pyna


class TestNSE:
    @pytest.fixture(scope="class")
    def pynet(self, reaclib_library):

        lib = reaclib_library.linking_nuclei(["p", "he4", "fe52",
                                              "co55", "ni56"])
        return pyna.NSENetwork(libraries=lib)

    @pytest.fixture(scope="class")
    def upper_net(self, tabular_library, reaclib_library):

        nuclei = ["p", "n", "cu72", "cu73", "ni72", "ni73"]
        tablib = tabular_library.linking_nuclei(nuclei)
        reaclib = reaclib_library.linking_nuclei(nuclei)

        lib = tablib + reaclib

        rates_to_remove = []
        for pair in lib.find_duplicate_links():
            for r in pair:
                if isinstance(r, pyna.rates.ReacLibRate):
                    rates_to_remove.append(r)
        for r in rates_to_remove:
            lib.remove_rate(r)

        return pyna.NSENetwork(libraries=lib)

    def test_nse_coul(self, pynet):

        rho = 1.0e7
        T = 6.0e9
        ye = 0.5

        nse_comp = pynet.get_comp_nse(rho, T, ye, use_coulomb_corr=True, use_unreliable_spins=False)

        nse_Xs = list(nse_comp.X.values())
        xsum = sum(nse_Xs)

        assert xsum == pytest.approx(1.0, rel=1.0e-10)

        expected = [
            0.009432340970989517, 0.4408337685107166, 0.007140861274572857,
            0.5187787534043365, 0.02381427583928699
        ]

        assert nse_Xs == pytest.approx(expected, rel=1.0e-10)

    def test_nse_no_coul(self, pynet):

        rho = 1.0e7
        T = 6.0e9
        ye = 0.5

        nse_comp = pynet.get_comp_nse(rho, T, ye, use_coulomb_corr=False, use_unreliable_spins=False)

        nse_Xs = list(nse_comp.X.values())
        xsum = sum(nse_Xs)

        assert xsum == pytest.approx(1.0, rel=1.0e-10)

        expected = [
            0.009096405176710216, 0.4630844653558956, 0.006705879462385591,
            0.5003022847190651, 0.020810965285946503
        ]

        assert nse_Xs == pytest.approx(expected, rel=1.0e-10)

    def test_nse_all_spin(self, upper_net):
        
        rho = 1.0e9
        T = 1.0e9
        ye = 0.45  # neutron-rich environment

        nse_comp = upper_net.get_comp_nse(rho, T, ye, use_coulomb_corr=True, use_unreliable_spins=True)

        nse_Xs = list(nse_comp.X.values())
        xsum = sum(nse_Xs)

        assert xsum == pytest.approx(1.0, rel=1.0e-10)

        assert nse_Xs == pytest.approx([0, 0.0790697674418601, 0, 0, 0.9209302325578685, 0], rel=1.0e-10)

    def test_nse_reliable_spin(self, upper_net):
        
        rho = 1.0e9
        T = 1.0e9
        ye = 0.45

        with pytest.raises(ValueError) as errorcode:
            upper_net.get_comp_nse(rho, T, ye, use_coulomb_corr=True, use_unreliable_spins=False)

        assert errorcode.value.args[0] == 'The spin of Ni73 is determined by a weak experimental' + \
            ' or theoretical argument. Pass in use_unreliable_spins=True as a parameter to override.'
