# This tests the PythonNetwork integrate_solution and NetworkSolution methods

import warnings

import numpy as np
import pytest
from numpy.testing import assert_allclose

from pynucastro import Composition, Nucleus, network_helper


class TestPythonIntegrate:
    @pytest.fixture(scope="class")
    @classmethod
    def net(cls):
        return network_helper(["he4", "c12"])

    @pytest.fixture(scope="class")
    @classmethod
    def sol(cls, net):
        rho = 1.e7
        T = 1.e9
        comp = Composition(net.unique_nuclei)
        comp.X[Nucleus("he4")] = 1.0
        tmax = 100.0

        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=UserWarning)
            sol = net.integrate_network(tmax, rho, T, comp.get_molar_array())

        return sol

    def test_sol_success(self, sol):
        assert sol.success

    def test_sol_X_at(self, sol):
        X = sol.X_at(100)

        assert_allclose(X, sol.X[:, -1], rtol=1.e-15, atol=1.e-30)
        assert_allclose(X, np.array([0.00216805, 0.99783195]), rtol=1.e-6, atol=1.e-30)

        Xs = sol.X_at([1, 10, 100])

        assert Xs.shape == (2, 3)

        assert_allclose(X, Xs[:, -1], rtol=1.e-15, atol=1.e-30)
