# unit tests for a rate collection
import pytest
from pytest import approx

import pynucastro as pyna


class TestRateCollectionJacobian:
    @pytest.fixture(scope="class")
    def rc(self, reaclib_library):
        mylib = reaclib_library.linking_nuclei(["he4", "c12", "o16", "ne20"],
                                               with_reverse=False)
        return pyna.RateCollection(libraries=[mylib])

    @pytest.fixture(scope="class")
    def comp(self, rc):
        _comp = pyna.Composition(rc.unique_nuclei)
        _comp.set_solar_like()
        return _comp

    def test_jac(self, rc, comp):
        rho = 1.e6
        T = 5.e8
        ymolar = comp.get_molar()

        jac = rc.evaluate_jacobian(rho, T, comp)

        # let's look now at jac(0, 0), which should be dYdot(He4)/dY(He4)
        #               and jac(1, 1), which should be dYdot(C12)/dY(C12)
        # there are 4 terms that go into this, which we will compute manually

        jac_alpha_alpha = 0.0
        jac_c12_c12 = 0.0

        # 3-alpha
        # r = rho**2 Y(alpha)**3 N_A<sigma v> / 6
        # 3 alphas are destroyed
        r = rc.get_rate_by_name("he4(aa,g)c12")
        jac_alpha_alpha += -3.0 * rho**2 * 3 * ymolar[pyna.Nucleus("he4")]**2 * r.eval(T) / 6.0
        # does not depend on c12

        # C12(a,g)O16
        # r = rho * Y(alpha) * Y(c12) N_A<sigma v>
        r = rc.get_rate_by_name("c12(a,g)o16")
        jac_alpha_alpha += -rho * ymolar[pyna.Nucleus("c12")] * r.eval(T)
        jac_c12_c12 += -rho * ymolar[pyna.Nucleus("he4")] * r.eval(T)

        # C12(C12,a)Ne20
        # r = rho * Y(c12)**2 N_A<sigma v> / 2
        r = rc.get_rate_by_name("c12(c12,a)ne20")
        # but the rate does not depend on alpha
        jac_c12_c12 += -2.0 * rho * 2.0 * ymolar[pyna.Nucleus("c12")] * r.eval(T) / 2.0

        # O16(a,g)Ne20
        # r = rho * Y(o16) * Y(alpha) N_A<sigma v>
        r = rc.get_rate_by_name("o16(a,g)ne20")
        jac_alpha_alpha += -rho * ymolar[pyna.Nucleus("o16")] * r.eval(T)
        # rate does not depend on c12

        assert jac_alpha_alpha == approx(jac[0, 0])
        assert jac_c12_c12 == approx(jac[1, 1])
