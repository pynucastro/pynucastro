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

        # scipy 1.15.0 rewrote MINPACK in C, which changed the results of
        # fsolve to roundoff. After compounding over multiple calls, the
        # difference is on the order of 1e-8 or so.
        expected_scipy_1_14 = [
            0.009432340970989385, 0.4408337685107169, 0.007140861274573644,
            0.5187787534044511, 0.023814275839292243
        ]
        expected_scipy_1_15 = [
            0.009432340953411016, 0.44083376971314875, 0.007140861228480336,
            0.5187787524376243, 0.02381427566734936
        ]
        assert (
            nse_Xs == pytest.approx(expected_scipy_1_14, rel=1.0e-10) or
            nse_Xs == pytest.approx(expected_scipy_1_15, rel=1.0e-10)
        )

    def test_nse_no_coul(self, pynet):

        rho = 1.0e7
        T = 6.0e9
        ye = 0.5

        nse_comp = pynet.get_comp_nse(rho, T, ye, use_coulomb_corr=False)

        nse_Xs = list(nse_comp.X.values())
        xsum = sum(nse_Xs)

        assert xsum == pytest.approx(1.0, rel=1.0e-10)

        expected_scipy_1_14 = [
            0.009096405176710995, 0.46308446535589587, 0.006705879462385593,
            0.5003022847190104, 0.02081096528594651
        ]
        expected_scipy_1_15 = [
            0.009096405158020583, 0.4630844666063867, 0.006705879416720939,
            0.5003022836911319, 0.02081096512773966
        ]
        assert (
            nse_Xs == pytest.approx(expected_scipy_1_14, rel=1.0e-10) or
            nse_Xs == pytest.approx(expected_scipy_1_15, rel=1.0e-10)
        )
