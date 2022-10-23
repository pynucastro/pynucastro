# unit tests for a rate collection
import pytest

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
        jac = rc.evaluate_jacobian(rho, T, comp)

        # let's look now at jac(0, 0), which should be dYdot(He4)/dY(He4)
        # there are 4 terms that go into this, which we will compute manually

        # 3-alpha

        # C12(a,g)O16

        # C12(C12,a)Ne20

        # O16(a,g)Ne20

