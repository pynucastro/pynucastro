"""Composition derivatives at fixed density and temperature."""

import numpy as np
import pytest

from pynucastro.constants import constants
from pynucastro.screening import (PlasmaState, ScreenFactors, chugunov_2007,
                                  chugunov_2009, debye_huckel,
                                  make_plasma_state_with_derivatives,
                                  potekhin_1998, screen5)

METHODS = [screen5, debye_huckel, chugunov_2007, chugunov_2009, potekhin_1998]

# set the composition properties for the plasma
YS = np.array([0.01, 0.05, 0.1, 0.025, 0.01])
ZS = np.array([0, 1, 2, 6, 8])


@pytest.mark.parametrize("method", METHODS)
@pytest.mark.parametrize("temp, dens", [(1.e9, 1.e3), (1.e8, 1.e5), (1.e8, 1.e6),
                                        (2.e8, 1.e9), (1.e6, 1.e5), (1.e5, 1.e12)])
def test_composition_derivatives(method, temp, dens):
    """Vary each Y_i independently, including a neutron and clipped regimes."""

    # the pair of nuclei don't directly affect the derivatives.  Here
    # we do He4 and C12.
    factors = ScreenFactors(2, 4, 6, 12)
    pstate = make_plasma_state_with_derivatives(temp, dens, YS, ZS)
    fscn_derivs = method(pstate, factors)

    # now compute the screening factor without the derivatives
    fscn = method(PlasmaState(temp, dens, YS, ZS), factors)

    assert fscn_derivs.value == pytest.approx(fscn, rel=2.e-13)

    # now compute the derivatives, ∂f/∂Yᵢ via the chain rule -- we can
    # get them for all species in one array expression.
    derivatives = fscn_derivs.deriv.d1 * ZS + fscn_derivs.deriv.d2 * ZS**2

    for i, abundance in enumerate(YS):
        # use a centered difference approximation to estimate the
        # derivatives
        step = 1.e-5 * abundance
        plus = YS.copy()
        minus = YS.copy()
        plus[i] += step
        minus[i] -= step
        expected = (method(PlasmaState(temp, dens, plus, ZS), factors) -
                    method(PlasmaState(temp, dens, minus, ZS), factors)) / (2 * step)
        assert derivatives[i] == pytest.approx(expected, rel=3.e-6, abs=2.e-6)

    if method in (chugunov_2007, potekhin_1998):
        assert fscn_derivs.deriv.d2 == 0.0


def test_state_derivatives():
    """Check analytic derivatives of precomputed plasma quantities."""
    temp, dens = 2.e8, 1.e9
    state = make_plasma_state_with_derivatives(temp, dens, YS, ZS)
    ordinary = PlasmaState(temp, dens, YS, ZS)

    # compute Y_e and S_2
    ye = np.sum(ZS * YS)
    s2 = np.sum(ZS**2 * YS)
    ytot = np.sum(YS)

    # these are the analytic derivatives of each component of
    # the plasma state with respect to Y_e and S_2
    expected = {"zbar": (1 / ytot, 0.0),
                "z2bar": (0.0, 1 / ytot),
                "n_e": (dens / constants.m_u_C18, 0.0),
                "aa": (ordinary.aa / (3 * ye), 0.0),
                "gamma_e_fac": (ordinary.gamma_e_fac / (3 * ye), 0.0),
                "qlam0z": (ordinary.qlam0z / (2 * (ye + s2)),) * 2}

    # make sure we get the expected values
    for name, (dye, ds2) in expected.items():
        field = getattr(state, name)
        assert field.value == pytest.approx(getattr(ordinary, name), rel=2.e-15)
        assert field.deriv.d1 == pytest.approx(dye)
        assert field.deriv.d2 == pytest.approx(ds2)
