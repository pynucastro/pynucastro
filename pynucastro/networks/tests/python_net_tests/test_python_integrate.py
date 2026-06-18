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

    @pytest.fixture(scope="class")
    @classmethod
    def sol_heating(cls, net):
        rho = 1.e7
        T = 1.e9
        comp = Composition(net.unique_nuclei)
        comp.X[Nucleus("he4")] = 1.0
        tmax = 100.0

        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=UserWarning)
            sol = net.integrate_network(tmax, rho, T, comp.get_molar_array(), self_heating=True)

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

    def test_sol_heating_X_at(self, sol_heating):

        X = sol_heating.X_at(100)

        assert_allclose(X, sol_heating.X[:, -1], rtol=1.e-15, atol=1.e-30)
        assert_allclose(X, np.array([0.36693869, 0.63306131]), rtol=1.e-6, atol=1.e-30)

        Xs = sol_heating.X_at([1, 10, 100])

        assert Xs.shape == (2, 3)

        assert_allclose(X, Xs[:, -1], rtol=1.e-15, atol=1.e-30)

    def test_sol_heating_T_at(self, sol_heating):

        T = sol_heating.T_at(100)

        assert T == pytest.approx(3.47848964e+09, rel=1.e-5, abs=1.e-10)
