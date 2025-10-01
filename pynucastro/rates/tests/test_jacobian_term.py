# unit tests for rates
import pytest
from pytest import approx

from pynucastro.networks import Composition
from pynucastro.nucdata import Nucleus


class TestJacTerm:

    @pytest.fixture(scope="class")
    def comp(self):
        nuc_list = [Nucleus("he4"), Nucleus("c12"), Nucleus("o16"),
                    Nucleus("ne20"), Nucleus("f20"), Nucleus("mg24")]
        comp = Composition(nuc_list)
        comp.set_solar_like()
        return comp

    def test_jac_3alpha(self, comp, reaclib_library):

        r = reaclib_library.get_rate_by_name("he4(aa,g)c12")

        rho = 1.e6
        T = 3.e8
        ymolar = comp.get_molar()

        # the full rate is Y(alpha)**3 rho**2 N_A <sigma v> / 6
        # where the 6 is the 3! that comes from like nuclei

        # consider drate/d(he4)

        # the derivative with respect to alpha is:
        # 3*Y(alpha)**2 rho**2 N_A <sigma v> / 6

        dr_dalpha = 3 * ymolar[Nucleus("he4")]**2 * rho**2 * r.eval(T) / 6
        assert r.eval_jacobian_term(T, rho, comp, Nucleus("he4")) == approx(dr_dalpha)

        # now consider the rate with respect to c12 -- the derivative
        # should be zero

        assert r.eval_jacobian_term(T, rho, comp, Nucleus("c12")) == 0.0

    def test_c12ag(self, comp, reaclib_library):

        r = reaclib_library.get_rate_by_name("c12(a,g)o16")

        rho = 1.e6
        T = 3.e8
        ymolar = comp.get_molar()

        # the full rate is Y(c12) Y(alpha) rho N_A <sigma v>

        # consider drate/d(he4)

        dr_dalpha = ymolar[Nucleus("c12")] * rho * r.eval(T)
        assert r.eval_jacobian_term(T, rho, comp, Nucleus("he4")) == approx(dr_dalpha)

        # now drate/d(c12)

        dr_dc12 = ymolar[Nucleus("he4")] * rho * r.eval(T)
        assert r.eval_jacobian_term(T, rho, comp, Nucleus("c12")) == approx(dr_dc12)

        # now drate/d(o16) should be 0
        assert r.eval_jacobian_term(T, rho, comp, Nucleus("o16")) == 0.0

    def test_tabular_rate(self, comp, tabular_library):

        r = tabular_library.get_rate_by_name("ne20(,)f20")

        rho = 5.e9
        T = 3.e8

        # this full rate is Y(ne20) lambda, where lambda is the 1/tau
        # read from the table

        # the rate does not dependent on alpha
        assert r.eval_jacobian_term(T, rho, comp, Nucleus("he4")) == 0.0

        # for dr/dY(ne20), we just have the raw rate from the table
        dr_dne20 = r.eval(T, rho=rho, comp=comp)
        assert r.eval_jacobian_term(T, rho, comp, Nucleus("ne20")) == dr_dne20
