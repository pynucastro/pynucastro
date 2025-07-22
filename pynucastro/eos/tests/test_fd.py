# unit tests for rates
from pytest import approx

from pynucastro.eos import FermiIntegral
from pynucastro.eos.difference_utils import fourth_order_diff

class TestFermiDirac:

    def test_k_0_5(self):
        # compare k = 1/2 solutions to the Gong et al. code

        k = 0.5

        eta = -100
        beta = 0.0

        f = FermiIntegral(k, eta, beta)
        f.evaluate()

        assert f.F == approx(3.2968314946796124e-044, abs=1.e-100, rel=1.e-15)
        assert f.dF_deta == approx(3.2968314946796139e-044, abs=1.e-100, rel=1.e-15)
        assert f.dF_dbeta == approx(1.2363118105048548e-044, abs=1.e-100, rel=1.e-15)
        assert f.d2F_deta2 == approx(3.2968314946796129e-044, abs=1.e-100, rel=1.e-15)
        assert f.d2F_detadbeta == approx(1.2363118105048556e-044, abs=1.e-100, rel=1.e-15)
        assert f.d2F_dbeta2 == approx(-7.7269488156553455E-045, abs=1.e-100, rel=1.e-15)

        eta = -50
        beta = 1.0

        f = FermiIntegral(k, eta, beta)
        f.evaluate()

        assert f.F == approx(2.2314386397062104E-022, abs=1.e-100, rel=1.e-15)
        assert f.dF_deta == approx(2.2314386397062108E-022, abs=1.e-100, rel=1.e-15)
        assert f.dF_dbeta == approx(4.4513616188460852E-023, abs=1.e-100, rel=1.e-15)
        assert f.d2F_deta2 == approx(2.2314386397062099E-022, abs=1.e-100, rel=1.e-15)
        assert f.d2F_detadbeta == approx(4.4513616188460870E-023, abs=1.e-100, rel=1.e-15)
        assert f.d2F_dbeta2 == approx(-1.0696566775877673E-023, abs=1.e-100, rel=1.e-15)

        eta = 1.0
        beta = 10.0

        f = FermiIntegral(k, eta, beta)
        f.evaluate()

        assert f.F == approx(4.3093121530612724, abs=1.e-100, rel=1.e-15)
        assert f.dF_deta == approx(3.0909917562756424, abs=1.e-100, rel=1.e-15)
        assert f.dF_dbeta == approx(0.19002958762134178, abs=1.e-100, rel=1.e-15)
        assert f.d2F_deta2 == approx(1.6796585258434078, abs=1.e-100, rel=1.e-15)
        assert f.d2F_detadbeta == approx(0.13977966489586524, abs=1.e-100, rel=1.e-15)
        assert f.d2F_dbeta2 == approx(-8.4816696742261236E-003, abs=1.e-100, rel=1.e-15)

        eta = 500.0
        beta = 100.0

        f = FermiIntegral(k, eta, beta)
        f.evaluate()

        assert f.F == approx(883930.45936891437, abs=1.e-100, rel=1.e-15)
        assert f.dF_deta == approx(3535.6046159037460, abs=1.e-100, rel=1.e-15)
        assert f.dF_dbeta == approx(4419.2988177917523, abs=1.e-100, rel=1.e-15)
        assert f.d2F_deta2 == approx(7.0710678132825819, abs=1.e-100, rel=1.e-15)
        assert f.d2F_detadbeta == approx(17.677315986879442, abs=1.e-100, rel=1.e-15)
        assert f.d2F_dbeta2 == approx(-22.094727366363866, abs=1.e-100, rel=1.e-15)

    def test_k_m0_5(self):
        # compare k = -1/2 solutions to the Gong et al. code

        k = -0.5

        eta = -100
        beta = 100.0

        f = FermiIntegral(k, eta, beta)
        f.evaluate()

        assert f.F == approx(2.7816742731063666e-043, abs=1.e-100, rel=1.e-15)
        assert f.dF_deta == approx(2.7816742731063678e-043, abs=1.e-100, rel=1.e-15)
        assert f.dF_dbeta == approx(1.2653970717385436e-045, abs=1.e-100, rel=1.e-15)
        assert f.d2F_deta2 == approx(2.7816742731063674e-043, abs=1.e-100, rel=1.e-15)
        assert f.d2F_detadbeta == approx(1.2653970717385444e-045, abs=1.e-100, rel=1.e-15)
        assert f.d2F_dbeta2 == approx(-5.9528644489672306e-048, abs=1.e-100, rel=1.e-15)

        eta = -50
        beta = 0.0

        f = FermiIntegral(k, eta, beta)
        f.evaluate()

        assert f.F == approx(3.4186200954570750e-022, abs=1.e-100, rel=1.e-15)
        assert f.dF_deta == approx(3.4186200954570741e-022, abs=1.e-100, rel=1.e-15)
        assert f.dF_dbeta == approx(4.2732751193213414e-023, abs=1.e-100, rel=1.e-15)
        assert f.d2F_deta2 == approx(3.4186200954570736e-022, abs=1.e-100, rel=1.e-15)
        assert f.d2F_detadbeta == approx(4.2732751193213450e-023, abs=1.e-100, rel=1.e-15)
        assert f.d2F_dbeta2 == approx(-1.6024781697455028e-023, abs=1.e-100, rel=1.e-15)

        eta = 100.0
        beta = 100.0

        f = FermiIntegral(k, eta, beta)
        f.evaluate()

        assert f.F == approx(707.87776608227455, abs=1.e-100, rel=1.e-15)
        assert f.dF_deta == approx(7.0717751162111329, abs=1.e-100, rel=1.e-14)
        assert f.dF_dbeta == approx(3.5323860528718423, abs=1.e-100, rel=1.e-15)
        assert f.d2F_deta2 == approx(-7.0773544184610238E-006, abs=1.e-100, rel=1.e-11)
        assert f.d2F_detadbeta == approx(3.5351802891431347E-002, abs=1.e-100, rel=1.e-15)
        assert f.d2F_dbeta2 == approx(-1.7633986737239839E-002, abs=1.e-100, rel=1.e-15)

    def test_k_1_5(self):
        # compare k = 3/2 solutions to the Gong et al. code

        k = 1.5

        eta = -75
        beta = 10.0

        f = FermiIntegral(k, eta, beta)
        f.evaluate()

        assert f.F == approx(1.2553979904636453E-032, abs=1.e-100, rel=1.e-15)
        assert f.dF_deta == approx(1.2553979904636447E-032, abs=1.e-100, rel=1.e-15)
        assert f.dF_dbeta == approx(5.7230338002235925E-034, abs=1.e-100, rel=1.e-15)
        assert f.d2F_deta2 == approx(1.2553979904636453E-032, abs=1.e-100, rel=1.e-15)
        assert f.d2F_detadbeta == approx(5.7230338002235968E-034, abs=1.e-100, rel=1.e-15)
        assert f.d2F_dbeta2 == approx(-2.6212332799199379E-035, abs=1.e-100, rel=1.e-15)

        eta = -20
        beta = 100.0

        f = FermiIntegral(k, eta, beta)
        f.evaluate()

        assert f.F == approx(2.9294159670904935E-008, abs=1.e-100, rel=1.e-15)
        assert f.dF_deta == approx(2.9294159663320433E-008, abs=1.e-100, rel=1.e-15)
        assert f.dF_dbeta == approx(1.4502712363335996E-010, abs=1.e-100, rel=1.e-15)
        assert f.d2F_deta2 == approx(2.9294159648151430E-008, abs=1.e-100, rel=1.e-15)
        assert f.d2F_detadbeta == approx(1.4502712359617485E-010, abs=1.e-100, rel=1.e-15)
        assert f.d2F_dbeta2 == approx(-7.1804917485691793E-013, abs=1.e-100, rel=1.e-15)

        eta = 40
        beta = 10000.0

        f = FermiIntegral(k, eta, beta)
        f.evaluate()

        assert f.F == approx(1517805.2872690351, abs=1.e-100, rel=1.e-15)
        assert f.dF_deta == approx(113369.99663886119, abs=1.e-100, rel=1.e-15)
        assert f.dF_dbeta == approx(75.889697517711070, abs=1.e-100, rel=1.e-15)
        assert f.d2F_deta2 == approx(5656.8613205601887, abs=1.e-100, rel=1.e-15)
        assert f.d2F_detadbeta == approx(5.6684715477425227, abs=1.e-100, rel=1.e-15)
        assert f.d2F_dbeta2 == approx(-3.7944565338813361E-003, abs=1.e-100, rel=1.e-15)

    def test_k_2_5(self):
        # compare k = 5/2 solutions to the Gong et al. code

        k = 2.5

        eta = -55
        beta = 15.0

        f = FermiIntegral(k, eta, beta)
        f.evaluate()

        assert f.F == approx(2.1821367305243183E-023, abs=1.e-100, rel=1.e-15)
        assert f.dF_deta == approx(2.1821367305243180E-023, abs=1.e-100, rel=1.e-15)
        assert f.dF_dbeta == approx(6.9671437808138560E-025, abs=1.e-100, rel=1.e-15)
        assert f.d2F_deta2 == approx(2.1821367305243186E-023, abs=1.e-100, rel=1.e-15)
        assert f.d2F_detadbeta == approx(6.9671437808138587E-025, abs=1.e-100, rel=1.e-15)
        assert f.d2F_dbeta2 == approx(-2.2261771808312711E-026, abs=1.e-100, rel=1.e-15)

        eta = 10
        beta = 1.e-4

        f = FermiIntegral(k, eta, beta)
        f.evaluate()

        assert f.F == approx(1034.9073854351607, abs=1.e-100, rel=1.e-15)
        assert f.dF_deta == approx(335.76592223107093, abs=1.e-100, rel=1.e-15)
        assert f.dF_dbeta == approx(2231.0524130570820, abs=1.e-100, rel=1.e-15)
        assert f.d2F_deta2 == approx(80.071134602962601, abs=1.e-100, rel=1.e-15)
        assert f.d2F_detadbeta == approx(905.09777777902468, abs=1.e-100, rel=1.e-15)
        assert f.d2F_dbeta2 == approx(-5201.1741381184684, abs=1.e-100, rel=1.e-15)

        eta = 100
        beta = 1.0

        f = FermiIntegral(k, eta, beta)
        f.evaluate()

        assert f.F == approx(17946771.869476113, abs=1.e-100, rel=1.e-15)
        assert f.dF_deta == approx(714843.05556091876, abs=1.e-100, rel=1.e-15)
        assert f.dF_dbeta == approx(8740888.8928929009, abs=1.e-100, rel=1.e-15)
        assert f.d2F_deta2 == approx(21361.250145668106, abs=1.e-100, rel=1.e-15)
        assert f.d2F_detadbeta == approx(350417.80107778130, abs=1.e-100, rel=1.e-15)
        assert f.d2F_dbeta2 == approx(-4257540.4021439729, abs=1.e-100, rel=1.e-15)

    def test_dfdeta(self):

        # test the derivative dF/dη computed via direct integration by
        # comparing against a difference approximation.

        eps = 1.e-8

        for k in [-0.5, 0.5, 1.5, 2.5]:
            for eta in [-70, 0, 50, 500]:
                for beta in [0, 30, 100]:
                    f0 = FermiIntegral(k, eta, beta)
                    f0.evaluate(do_second_derivs=False)

                    eta_new = eta * (1.0 + eps)
                    if eta == 0.0:
                        eta_new = eps

                    f1 = FermiIntegral(k, eta_new, beta)
                    f1.evaluate(do_first_derivs=False, do_second_derivs=False)

                    deriv = (f1.F - f0.F) / (eta_new - eta)

                    assert f0.dF_deta == approx(deriv, abs=1.e-100, rel=1.e-6)

    def test_dfdbeta(self):

        # test the derivative dF/dβ computed via direct integration by
        # comparing against a difference approximation.

        eps = 1.e-8

        for k in [-0.5, 0.5, 1.5, 2.5]:
            for eta in [-70, 0, 50, 500]:
                for beta in [0, 30, 100]:
                    f0 = FermiIntegral(k, eta, beta)
                    f0.evaluate(do_second_derivs=False)

                    beta_new = beta * (1.0 + eps)
                    if beta == 0.0:
                        beta_new = eps

                    f1 = FermiIntegral(k, eta, beta_new)
                    f1.evaluate(do_first_derivs=False, do_second_derivs=False)

                    deriv = (f1.F - f0.F) / (beta_new - beta)

                    assert f0.dF_dbeta == approx(deriv, abs=1.e-100, rel=1.e-6)

    def test_d2fdeta2(self):

        # test the derivative d^2F/dη^2 computed via direct
        # integration by comparing against a difference approximation.

        eps = 1.e-8

        for k in [-0.5, 0.5, 1.5, 2.5]:
            for eta in [-50, 0, 250]:
                for beta in [1, 25, 200]:
                    f0 = FermiIntegral(k, eta, beta)
                    f0.evaluate()

                    deta = abs(eps * eta)
                    if deta == 0.0:
                        deta = eps

                    def _kernel(_eta):
                        _f = FermiIntegral(k, _eta, beta)
                        _f.evaluate(do_second_derivs=False)
                        return _f.dF_deta

                    deriv = fourth_order_diff(_kernel, eta, deta)

                    assert f0.d2F_deta2 == approx(deriv, abs=1.e-100, rel=5.e-4)

    def test_d2fdbeta2(self):

        # test the derivative d^2F/dβ^2 computed via direct
        # integration by comparing against a difference approximation.

        eps = 1.e-8

        for k in [-0.5, 0.5, 1.5, 2.5]:
            for eta in [-50, 0, 250]:
                for beta in [1, 25, 200]:
                    f0 = FermiIntegral(k, eta, beta)
                    f0.evaluate()

                    dbeta = abs(eps * beta)
                    if dbeta == 0.0:
                        dbeta = eps

                    def _kernel(_beta):
                        _f = FermiIntegral(k, eta, _beta)
                        _f.evaluate(do_second_derivs=False)
                        return _f.dF_bdeta

                    deriv = fourth_order_diff(_kernel, beta, dbeta)

                    assert f0.d2F_dbeta2 == approx(deriv, abs=1.e-100, rel=1.e-6)
