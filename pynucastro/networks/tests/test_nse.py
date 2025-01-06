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
        assert nse_Xs[0] == pytest.approx(0.00943234095341063, rel=1.0e-10)
        assert nse_Xs[1] == pytest.approx(0.44083376851072315, rel=1.0e-10)
        assert nse_Xs[2] == pytest.approx(0.007140861274574405, rel=1.0e-10)
        assert nse_Xs[3] == pytest.approx(0.5187787534043939, rel=1.0e-10)
        assert nse_Xs[4] == pytest.approx(0.023814275839289616, rel=1.0e-10)

    def test_nse_no_coul(self, pynet):

        rho = 1.0e7
        T = 6.0e9
        ye = 0.5

        nse_comp = pynet.get_comp_nse(rho, T, ye, use_coulomb_corr=False)

        nse_Xs = list(nse_comp.X.values())
        xsum = sum(nse_Xs)

        assert xsum == pytest.approx(1.0, rel=1.0e-10)
        assert nse_Xs[0] == pytest.approx(0.009096405158020583, rel=1.0e-10)
        assert nse_Xs[1] == pytest.approx(0.46308446535589587, rel=1.0e-10)
        assert nse_Xs[2] == pytest.approx(0.006705879462385593, rel=1.0e-10)
        assert nse_Xs[3] == pytest.approx(0.5003022847190655, rel=1.0e-10)
        assert nse_Xs[4] == pytest.approx(0.02081096528594651, rel=1.0e-10)
