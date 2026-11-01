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
    def test_jacobian_no_screening(self, net, temp, dens):

        comp=Composition(net.unique_nuclei, init="uniform")
        state = ThermoState(rho=dens, T=temp, comp=comp)

        jac_analytic = net.evaluate_jacobian(state)

        jac_numerical = np.zeros_like(jac_analytic)

        # do a finite-difference.  First compute the unperturbed
        # Ydots
        dydt0 = net.evaluate_ydots(state)

        h = 1.e-8

        for jcol, jnuc in enumerate(net.unique_nuclei):

            # perturb the current nucleus
            cnuc_old = comp.X[jnuc]
            comp.X[jnuc] += h

            state_pert = ThermoState(rho=dens, T=temp, comp=comp)
            dydt_pert = net.evaluate_ydots(state_pert)

            comp.X[jnuc] = cnuc_old

            for irow, inuc in enumerate(net.unique_nuclei):
                # we perturbed X, so we need to convert d/dX to d/dY
                jac_numerical[irow, jcol] = (dydt_pert[inuc] - dydt0[inuc]) / h * jnuc.A

        # compute the relative error
        err = np.abs(jac_numerical - jac_analytic)

        rerr = np.divide(err, jac_analytic, where=(jac_analytic != 0), out=None)

        assert rerr.max() < 1.e-2

