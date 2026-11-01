"""Test the RateCollection Jacobian against a finite-difference
approximation.

This is a simple net that has some weak rates and some approximate
rates, so it should test the composition derivative embedded in the
approximation and the composition dependence via Ye from the weak
rates.

"""

import numpy as np
import pytest

from pynucastro.networks import RateCollection
from pynucastro.nucdata import Composition
from pynucastro.rates import DerivedRate, ThermoState
from pynucastro.screening import chugunov_2007


class TestJacobianFD:
    @pytest.fixture(scope="class")
    @classmethod
    def net(cls, reaclib_library, langanke_library):

        nuc = ["p", "n", "he4",
               "mn51",
               "fe52", "fe53", "fe54", "fe55", "fe56",
               "co55", "co56", "co57",
               "ni56", "ni57", "ni58", "cu59", "zn60"]

        lib = reaclib_library.linking_nuclei(nuc)
        tlib = langanke_library.linking_nuclei(nuc)
        lib += tlib

        lib.eliminate_duplicates()

        for r in lib.backward().get_rates():
            if fr := lib.get_rate_by_nuclei(r.products, r.reactants):
                lib.remove_rate(r)
                lib.add_rate(DerivedRate(fr, use_pf=True))

        net = RateCollection(libraries=[lib])

        net.make_nn_g_approx(intermediate_nuclei=["fe53", "fe55", "ni57"])
        net.remove_nuclei(["fe53", "fe55", "ni57"])

        return net

    @pytest.mark.parametrize("temp, dens", [(1.e8, 1.e5), (1.e9, 1.e6), (4.e9, 1.e7)])
    @pytest.mark.parametrize("screen_method", [None, chugunov_2007], ids=["no_screening", "chugunov_2007"])
    def test_jacobian(self, net, temp, dens, screen_method):

        comp = Composition(net.unique_nuclei, init="uniform")
        state = ThermoState(rho=dens, T=temp, comp=comp)

        jac_analytic = net.evaluate_jacobian(state, screen_func=screen_method)

        jac_numerical = np.zeros_like(jac_analytic)

        # Estimate roundoff from all production and destruction terms,
        # since their sum can be much larger than the net Ydot.
        fluxes = net.evaluate_rates(state, screen_func=screen_method)
        flux_scale = np.array([
            sum(r.reactant_count(n) * abs(fluxes[r]) for r in net.nuclei_consumed[n]) +
            sum(r.product_count(n) * abs(fluxes[r]) for r in net.nuclei_produced[n])
            for n in net.unique_nuclei
        ])

        h = 1.e-5

        for jcol, jnuc in enumerate(net.unique_nuclei):

            # perturb the current nucleus
            cnuc_old = comp.X[jnuc]
            comp.X[jnuc] = cnuc_old + h

            state_pert = ThermoState(rho=dens, T=temp, comp=comp)
            dydt_plus = net.evaluate_ydots(state_pert, screen_func=screen_method)

            comp.X[jnuc] = cnuc_old - h
            dydt_minus = net.evaluate_ydots(state_pert, screen_func=screen_method)

            comp.X[jnuc] = cnuc_old

            for irow, inuc in enumerate(net.unique_nuclei):
                # we perturbed X, so we need to convert d/dX to d/dY
                jac_numerical[irow, jcol] = (dydt_plus[inuc] - dydt_minus[inuc]) / (2 * h) * jnuc.A

        # Allow 1% relative error plus an estimate of the roundoff in the
        # finite difference.  Each Ydot evaluation has rounding errors on
        # the scale of eps * flux_scale, where eps is machine precision.
        # Use the sum of absolute production and destruction terms, since
        # cancellation can make the net Ydot much smaller than these terms.
        # The two evaluations contribute roughly 2 * eps * flux_scale of
        # uncertainty to the numerator; dividing by 2*h leaves a scale of
        # eps * flux_scale / h.  Thus a smaller h amplifies roundoff.
        # Multiplying by A_j converts this d/dX_j uncertainty to d/dY_j,
        # giving the same inverse-second units as the Jacobian entries.
        # The factor of 10 is a practical allowance for flux evaluation and
        # summation errors, not a rigorous bound.  Entries below this floor
        # cannot be reliably checked to 1% by this finite difference.
        err = np.abs(jac_numerical - jac_analytic)
        mass_numbers = np.array([n.A for n in net.unique_nuclei])
        roundoff_floor = 10 * np.finfo(float).eps * flux_scale[:, None] * mass_numbers[None, :] / h
        error_bound = 1.e-2 * np.abs(jac_analytic) + roundoff_floor

        # argmax returns the flat index of the entry with the largest
        # excess error over its bound.  unravel_index converts that index
        # into a row and column so the failure message can name the nuclei.
        irow, jcol = np.unravel_index(np.argmax(err - error_bound), err.shape)
        assert np.all(err <= error_bound), (
            f"Jacobian entry ({net.unique_nuclei[irow]}, {net.unique_nuclei[jcol]}): "
            f"analytic={jac_analytic[irow, jcol]}, numerical={jac_numerical[irow, jcol]}, "
            f"error={err[irow, jcol]}, bound={error_bound[irow, jcol]}"
        )
