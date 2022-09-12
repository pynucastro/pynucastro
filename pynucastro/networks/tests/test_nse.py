import pynucastro as pyna
import pytest


class TestNSE:
    @pytest.fixture(scope="class")
    def pynet(self, reaclib_library):

        lib = reaclib_library.linking_nuclei(["p", "he4", "fe52",
                                              "co55", "ni56"])
        return pyna.RateCollection(libraries=lib)

    def test_nse_coul(self, pynet):

        rho = 1.0e7
        T = 6.0e9
        ye = 0.5

        nse_comp = pynet.get_comp_NSE(rho, T, ye, use_coulomb_corr=True)
        nse_Xs = list(nse_comp.X.values())

        assert nse_Xs[0] == pytest.approx(0.00573778851016, rel=1.0e-10)
        assert nse_Xs[1] == pytest.approx(0.50055406187850, rel=1.0e-10)
        assert nse_Xs[2] == pytest.approx(0.03711671370420, rel=1.0e-10)
        assert nse_Xs[3] == pytest.approx(0.31557836805821, rel=1.0e-10)
        assert nse_Xs[4] == pytest.approx(0.14101306784484, rel=1.0e-10)

    def test_nse_no_coul(self, pynet):

        rho = 1.0e7
        T = 6.0e9
        ye = 0.5

        nse_comp = pynet.get_comp_NSE(rho, T, ye, use_coulomb_corr=False)
        nse_Xs = list(nse_comp.X.values())

        assert nse_Xs[0] == pytest.approx(0.00558282578105, rel=1.0e-10)
        assert nse_Xs[1] == pytest.approx(0.52648694688751, rel=1.0e-10)
        assert nse_Xs[2] == pytest.approx(0.03543576484397, rel=1.0e-10)
        assert nse_Xs[3] == pytest.approx(0.30705541795776, rel=1.0e-10)
        assert nse_Xs[4] == pytest.approx(0.12543904452967, rel=1.0e-10)
