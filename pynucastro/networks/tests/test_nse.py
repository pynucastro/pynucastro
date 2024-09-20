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
        assert nse_Xs[0] == pytest.approx(0.009432537966997221, rel=1.0e-10)
        assert nse_Xs[1] == pytest.approx(0.44084124187076745, rel=1.0e-10)
        assert nse_Xs[2] == pytest.approx(0.007118462338888465, rel=1.0e-10)
        assert nse_Xs[3] == pytest.approx(0.5187895881848553, rel=1.0e-10)
        assert nse_Xs[4] == pytest.approx(0.02381816963849376, rel=1.0e-10)

    def test_nse_no_coul(self, pynet):

        rho = 1.0e7
        T = 6.0e9
        ye = 0.5

        nse_comp = pynet.get_comp_nse(rho, T, ye, use_coulomb_corr=False)

        nse_Xs = list(nse_comp.X.values())
        xsum = sum(nse_Xs)

        assert xsum == pytest.approx(1.0, rel=1.0e-10)
        assert nse_Xs[0] == pytest.approx(0.009096582765040962, rel=1.0e-10)
        assert nse_Xs[1] == pytest.approx(0.46309222580780524, rel=1.0e-10)
        assert nse_Xs[2] == pytest.approx(0.006684828034485319, rel=1.0e-10)
        assert nse_Xs[3] == pytest.approx(0.5003120520772527, rel=1.0e-10)
        assert nse_Xs[4] == pytest.approx(0.020814311315418953, rel=1.0e-10)
