import pytest

import pynucastro as pyna


class TestNSE:
    @pytest.fixture(scope="class")
    def pynet(self, reaclib_library):

        lib = reaclib_library.linking_nuclei(["p", "he4", "fe52",
                                              "co55", "ni56"])
        return pyna.NSENetwork(libraries=lib)

    def test_nse_coul(self, pynet):

        rho = 1.0e7
        T = 6.0e9
        ye = 0.5

        nse_comp = pynet.get_comp_nse(rho, T, ye, use_coulomb_corr=True)

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

        nse_comp = pynet.get_comp_nse(rho, T, ye, use_coulomb_corr=False)

        nse_Xs = list(nse_comp.X.values())
        xsum = sum(nse_Xs)

        assert xsum == pytest.approx(1.0, rel=1.0e-10)

        expected = [
            0.009096405176710216, 0.4630844653558956, 0.006705879462385591,
            0.5003022847190651, 0.020810965285946503
        ]

        assert nse_Xs == pytest.approx(expected, rel=1.0e-10)
