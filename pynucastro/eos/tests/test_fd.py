# unit tests for rates
import math

import pytest
from pytest import approx

from pynucastro.eos import FermiIntegrals

class TestFermiDirac:

    def test_k_0_5(self):
        # compare k = 1/2 solutions to the Gong et al. code

        k = 0.5

        eta = -100
        beta = 0.0

        f = FermiIntegrals(k, eta, beta)
        f.evaluate()

        assert f.F == approx(3.2968314946796124e-044, abs=1.e-100, rel=1.e-13)
        assert f.dF_deta == approx(3.2968314946796139e-044, abs=1.e-100, rel=1.e-13)
        assert f.dF_dbeta == approx(1.2363118105048548e-044, abs=1.e-100, rel=1.e-13)
        assert f.d2F_deta2 == approx(3.2968314946796129e-044, abs=1.e-100, rel=1.e-13)
        assert f.d2F_detadbeta == approx(1.2363118105048556e-044, abs=1.e-100, rel=1.e-13)
        assert f.d2F_dbeta2 == approx(-7.7269488156553455E-045, abs=1.e-100, rel=1.e-13)

        eta = -50
        beta = 1.0

        f = FermiIntegrals(k, eta, beta)
        f.evaluate()

        assert f.F == approx(2.2314386397062104E-022, abs=1.e-100, rel=1.e-13)
        assert f.dF_deta == approx(2.2314386397062108E-022, abs=1.e-100, rel=1.e-13)
        assert f.dF_dbeta == approx(4.4513616188460852E-023, abs=1.e-100, rel=1.e-13)
        assert f.d2F_deta2 == approx(2.2314386397062099E-022, abs=1.e-100, rel=1.e-13)
        assert f.d2F_detadbeta == approx(4.4513616188460870E-023, abs=1.e-100, rel=1.e-13)
        assert f.d2F_dbeta2 == approx(-1.0696566775877673E-023, abs=1.e-100, rel=1.e-13)

        eta = 1.0
        beta = 10.0

        f = FermiIntegrals(k, eta, beta)
        f.evaluate()

        assert f.F == approx(4.3093121530612724, abs=1.e-100, rel=1.e-13)
        assert f.dF_deta == approx(3.0909917562756424, abs=1.e-100, rel=1.e-13)
        assert f.dF_dbeta == approx(0.19002958762134178, abs=1.e-100, rel=1.e-13)
        assert f.d2F_deta2 == approx(1.6796585258434078, abs=1.e-100, rel=1.e-13)
        assert f.d2F_detadbeta == approx(0.13977966489586524, abs=1.e-100, rel=1.e-13)
        assert f.d2F_dbeta2 == approx(-8.4816696742261236E-003, abs=1.e-100, rel=1.e-13)

