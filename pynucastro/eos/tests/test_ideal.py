# test the ideal gas

from pytest import approx

from pynucastro.constants import constants
from pynucastro.eos import IdealGasEOS
from pynucastro.nucdata import Composition, Nucleus


class TestIdealGasEOS:

    def test_ionized(self):

        # test fully-ionized He4
        comp = Composition([Nucleus("he4")])
        comp.X[Nucleus("he4")] = 1.0

        # the total mean molecular weight for ionized He4 is 4/3
        mu = 4.0 / 3.0

        rho, T = 1.e4, 1.e7
        state = IdealGasEOS(include_electrons=True).pe_state(rho, T, comp)

        # compute the pressure manually
        p = rho * constants.k * T / (mu * constants.m_u)

        assert state.p == approx(p, rel=1.e-12)
        assert state.e == approx(1.5 * p / rho, rel=1.e-12)
        assert state.dp_drho == approx(p / rho, rel=1.e-12)
