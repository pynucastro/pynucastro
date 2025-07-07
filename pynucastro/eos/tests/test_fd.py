# unit tests for rates
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

        eta = 500.0
        beta = 100.0

        f = FermiIntegrals(k, eta, beta)
        f.evaluate()

        assert f.F == approx(883930.45936891437, abs=1.e-100, rel=1.e-13)
        assert f.dF_deta == approx(3535.6046159037460, abs=1.e-100, rel=1.e-13)
        assert f.dF_dbeta == approx(4419.2988177917523, abs=1.e-100, rel=1.e-13)
        assert f.d2F_deta2 == approx(7.0710678132825819, abs=1.e-100, rel=1.e-13)
        assert f.d2F_detadbeta == approx(17.677315986879442, abs=1.e-100, rel=1.e-13)
        assert f.d2F_dbeta2 == approx(-22.094727366363866, abs=1.e-100, rel=1.e-13)

    def test_k_m0_5(self):
        # compare k = -1/2 solutions to the Gong et al. code

        k = -0.5

        eta = -100
        beta = 100.0

        f = FermiIntegrals(k, eta, beta)
        f.evaluate()

        assert f.F == approx(2.7816742731063666e-043, abs=1.e-100, rel=1.e-13)
        assert f.dF_deta == approx(2.7816742731063678e-043, abs=1.e-100, rel=1.e-13)
        assert f.dF_dbeta == approx(1.2653970717385436e-045, abs=1.e-100, rel=1.e-13)
        assert f.d2F_deta2 == approx(2.7816742731063674e-043, abs=1.e-100, rel=1.e-13)
        assert f.d2F_detadbeta == approx(1.2653970717385444e-045, abs=1.e-100, rel=1.e-13)
        assert f.d2F_dbeta2 == approx(-5.9528644489672306e-048, abs=1.e-100, rel=1.e-13)

        eta = -50
        beta = 0.0

        f = FermiIntegrals(k, eta, beta)
        f.evaluate()

        assert f.F == approx(3.4186200954570750e-022, abs=1.e-100, rel=1.e-13)
        assert f.dF_deta == approx(3.4186200954570741e-022, abs=1.e-100, rel=1.e-13)
        assert f.dF_dbeta == approx(4.2732751193213414e-023, abs=1.e-100, rel=1.e-13)
        assert f.d2F_deta2 == approx(3.4186200954570736e-022, abs=1.e-100, rel=1.e-13)
        assert f.d2F_detadbeta == approx(4.2732751193213450e-023, abs=1.e-100, rel=1.e-13)
        assert f.d2F_dbeta2 == approx(-1.6024781697455028e-023, abs=1.e-100, rel=1.e-13)

        eta = 100.0
        beta = 100.0

        f = FermiIntegrals(k, eta, beta)
        f.evaluate()

        assert f.F == approx(707.87776608227455, abs=1.e-100, rel=1.e-13)
        assert f.dF_deta == approx(7.0717751162111329, abs=1.e-100, rel=1.e-13)
        assert f.dF_dbeta == approx(3.5323860528718423, abs=1.e-100, rel=1.e-13)
        assert f.d2F_deta2 == approx(-7.0773544184610238E-006, abs=1.e-100, rel=1.e-13)
        assert f.d2F_detadbeta == approx(3.5351802891431347E-002, abs=1.e-100, rel=1.e-13)
        assert f.d2F_dbeta2 == approx(-1.7633986737239839E-002, abs=1.e-100, rel=1.e-13)

    def test_dfdeta(self):

        # test the derivative dF/dβ computed via direct integration by
        # comparing against a difference approximation.

        eps = 1.e-8

        for k in [-0.5, 0.5, 1.5, 2.5]:
            for eta in [-70, 0, 50, 500]:
                for beta in [0, 30, 100]:
                    f0 = FermiIntegrals(k, eta, beta)
                    f0.evaluate()

                    eta_new = eta * (1.0 + eps)
                    if eta == 0.0:
                        eta_new = eps

                    f1 = FermiIntegrals(k, eta_new, beta)
                    f1.evaluate()

                    deriv = (f1.F - f0.F) / (eta_new - eta)

                    assert f0.dF_deta == approx(deriv, abs=1.e-100, rel=1.e-6)

    def test_d2fdeta2(self):

        # test the derivative d^2F/dβ^2 computed via direct
        # integration by comparing against a difference approximation.

        eps = 1.e-8

        for k in [-0.5, 0.5, 1.5, 2.5]:
            for eta in [-50, 0, 250]:
                for beta in [1, 25, 200]:
                    f0 = FermiIntegrals(k, eta, beta)
                    f0.evaluate()

                    deta = abs(eps * eta)
                    if deta == 0.0:
                        deta = eps

                    fp = FermiIntegrals(k, eta+deta, beta)
                    fp.evaluate()

                    fp2 = FermiIntegrals(k, eta+2*deta, beta)
                    fp2.evaluate()

                    fm = FermiIntegrals(k, eta-deta, beta)
                    fm.evaluate()

                    fm2 = FermiIntegrals(k, eta-2*deta, beta)
                    fm2.evaluate()

                    #deriv2 = (fp.F - 2.0*f0.F + fm.F) / deta**2
                    #deriv2 = (-fm2.F + 16*fm.F - 30*f0.F + 16*fp.F - fp2.F) / 12 / deta**2
                    deriv2 = 0.5 * (fp.dF_deta - fm.dF_deta) / deta

                    assert f0.d2F_deta2 == approx(deriv2, abs=1.e-100, rel=1.e-3)
