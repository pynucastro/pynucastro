import pytest

import pynucastro as pyna


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

        nse_comp = pynet.get_comp_nse(rho, T, ye, use_coulomb_corr=True)

        sum = 0.0
        for nuc in nse_comp.X:
            sum += nse_comp.X[nuc]

        nse_Xs = list(nse_comp.X.values())

        assert sum == pytest.approx(1.0, rel=1.0e-10)
        assert nse_Xs[0] == pytest.approx(0.009432528259141131, rel=1.0e-10)
        assert nse_Xs[1] == pytest.approx(0.44084190707382365, rel=1.0e-10)
        assert nse_Xs[2] == pytest.approx(0.007118438690709681, rel=1.0e-10)
        assert nse_Xs[3] == pytest.approx(0.5187890542527211, rel=1.0e-10)
        assert nse_Xs[4] == pytest.approx(0.023818071723544595, rel=1.0e-10)

    def test_nse_no_coul(self, pynet):

        rho = 1.0e7
        T = 6.0e9
        ye = 0.5

        nse_comp = pynet.get_comp_nse(rho, T, ye, use_coulomb_corr=False)

        sum = 0.0
        for nuc in nse_comp.X:
            sum += nse_comp.X[nuc]

        nse_Xs = list(nse_comp.X.values())

        assert sum == pytest.approx(1.0, rel=1.0e-10)
        assert nse_Xs[0] == pytest.approx(0.009096572438869752, rel=1.0e-10)
        assert nse_Xs[1] == pytest.approx(0.4630929175748844, rel=1.0e-10)
        assert nse_Xs[2] == pytest.approx(0.006684804511147247, rel=1.0e-10)
        assert nse_Xs[3] == pytest.approx(0.5003114841378136, rel=1.0e-10)
        assert nse_Xs[4] == pytest.approx(0.020814221337257847, rel=1.0e-10)
