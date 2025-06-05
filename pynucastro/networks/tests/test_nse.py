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
            0.0094323409692503, 0.4408337686296852, 0.007140861270012667,
            0.518778753308708, 0.023814275822275126
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
            0.00909640517486066, 0.46308446547962107, 0.0067058794578672525,
            0.5003022846173812, 0.020810965270293028
        ]

        assert nse_Xs == pytest.approx(expected, rel=1.0e-10)
