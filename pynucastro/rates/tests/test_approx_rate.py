# unit tests for rates

import pytest
from pytest import approx

import pynucastro as pyna
from pynucastro.rates import create_double_neutron_capture


class TestAlphaGammaTfactors:
    @pytest.fixture(scope="class")
    def rc(self, reaclib_library):
        mylib = reaclib_library.linking_nuclei(["mg24", "al27", "si28", "he4", "p"])
        return pyna.RateCollection(libraries=[mylib])

    @pytest.fixture(scope="class")
    def rp(self, rc):
        return rc.get_rate("he4_mg24_to_si28")

    @pytest.fixture(scope="class")
    def rs(self, rc):
        return [rc.get_rate("he4_mg24_to_p_al27"), rc.get_rate("p_al27_to_si28")]

    @pytest.fixture(scope="class")
    def ar(self, rc, rp, rs):
        rp_reverse = rc.get_rate("si28_to_he4_mg24")
        rs_reverse = [rc.get_rate("si28_to_p_al27"), rc.get_rate("p_al27_to_he4_mg24")]

        # approximate Mg24(a,g)Si28 together with Mg24(a,p)Al27(p,g)Si28

        return pyna.ApproximateRate(primary_rate=rp, secondary_rates=rs,
                                    primary_reverse=rp_reverse, secondary_reverse=rs_reverse)

    def test_label(self, ar):
        assert ar.fname == "Mg24_He4_to_Si28_approx"

    def test_low_temp(self, ar, rp, rs):
        # at low temperatures, the approximate (a,g) should be ~ (a,g) + (a,p)

        T = 5.e8
        assert ar.eval(T) == approx(rp.eval(T) + rs[0].eval(T), 1.e-6)

        T = 1.e9
        assert ar.eval(T) == approx(rp.eval(T) + rs[0].eval(T), 0.1)

    def test_child_rates(self, ar):

        cr = ar.get_child_rates()
        assert cr[0].fname == "He4_Mg24_to_Si28"
        assert cr[1].fname == "He4_Mg24_to_p_Al27"
        assert cr[2].fname == "p_Al27_to_Si28"
        assert cr[3].fname == "Si28_to_He4_Mg24"
        assert cr[4].fname == "Si28_to_p_Al27"
        assert cr[5].fname == "p_Al27_to_He4_Mg24"

        assert len(cr) == 6


class TestDoubleN:

    @pytest.fixture(scope="class")
    def nn(self, reaclib_library):
        return create_double_neutron_capture(reaclib_library, "Fe52", "Fe54")

    def test_nn_name(self, nn):
        """ this is run before each test """

        rf, rr = nn

        assert rf.fname == "Fe52_n_n_to_Fe54_approx"
        assert rr.fname == "Fe54_to_Fe52_n_n_approx"

    def test_Q(self, nn):

        rf, _ = nn

        assert rf.Q == approx(24.06513619999896)

    def test_func(self, nn):

        rf, _ = nn

        assert rf.function_string_py() == \
"""@numba.njit()
def Fe52_n_n_to_Fe54_approx(rate_eval, tf, rho=None, Y=None):
    Yn = Y[jn]
    r1_ng = rate_eval.n_Fe52_to_Fe53
    r2_ng = rate_eval.n_Fe53_to_Fe54
    r1_gn = rate_eval.Fe53_to_n_Fe52
    rate = r1_ng * r2_ng / (rho * Yn * r2_ng + r1_gn)
    rate_eval.Fe52_n_n_to_Fe54_approx = rate

"""

    def test_particle_factor(self, nn):

        rf, _ = nn

        assert not rf.use_identical_particle_factor

    def test_eval(self, nn):

        rf, rr = nn

        comp = pyna.Composition(["n", "fe52", "fe54"])
        comp.set_equal()

        T = 1.e9
        rho = 2.e7

        rf1, _, __, ___ = rr.hidden_rates

        assert rf1.eval(T, rho=rho, comp=comp) == approx(6511141.861787519)

        Yn = comp.get_molar()[pyna.Nucleus("n")]

        # at low temperature, the reverse rate is ~ zero,
        # so the total rate will just be the first forward
        # rate scaled by rho Y(n)

        assert rf.eval(T, rho=rho, comp=comp) == approx(rf1.eval(T, rho=rho, comp=comp) / rho / Yn)
