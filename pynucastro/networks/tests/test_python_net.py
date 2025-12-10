# unit tests for rates
import pytest

from pynucastro import rates


class TestPythonNetwork:
    @pytest.fixture(scope="class")
    def rate1(self):
        return rates.ReacLibRate("c13-pg-n14-nacr")

    @pytest.fixture(scope="class")
    def rate2(self):
        return rates.ReacLibRate("he4-pphe3-he3-nacr")

    @pytest.fixture(scope="class")
    def rate3(self):
        return rates.ReacLibRate("he4-npahe3-li7-mafo")

    @pytest.fixture(scope="class")
    def rate4(self):
        return rates.ReacLibRate("p-p-d-ec")

    def test_ydot_string(self, rate1, rate2, rate3, rate4):
        ydot1 = rate1.ydot_string_py()
        assert ydot1 == "rho*Y[jp]*Y[jc13]*rate_eval.p_C13_to_N14"

        ydot2 = rate2.ydot_string_py()
        assert ydot2 == "5.00000000000000e-01*rho**2*Y[jp]**2*Y[jhe4]*rate_eval.p_p_He4_to_He3_He3"

        ydot3 = rate3.ydot_string_py()
        assert ydot3 == "5.00000000000000e-01*rho**3*Y[jn]*Y[jp]*Y[jhe4]**2*rate_eval.n_p_He4_He4_to_He3_Li7"

        ydot4 = rate4.ydot_string_py()
        assert ydot4 == "5.00000000000000e-01*rho**2*ye(Y)*Y[jp]**2*rate_eval.p_p_to_d_weak_electron_capture"

    def test_jacobian_string(self, rate1, rate2, rate3, rate4):
        jac1 = rate1.jacobian_string_py(rate1.reactants[0])
        assert jac1 == "rho*Y[jc13]*rate_eval.p_C13_to_N14"

        jac2 = rate2.jacobian_string_py(rate2.reactants[0])
        assert jac2 == "5.00000000000000e-01*rho**2*2*Y[jp]*Y[jhe4]*rate_eval.p_p_He4_to_He3_He3"

        jac3 = rate3.jacobian_string_py(rate3.reactants[0])
        assert jac3 == "5.00000000000000e-01*rho**3*Y[jp]*Y[jhe4]**2*rate_eval.n_p_He4_He4_to_He3_Li7"

        jac4 = rate4.jacobian_string_py(rate4.reactants[0])
        assert jac4 == "5.00000000000000e-01*rho**2*ye(Y)*2*Y[jp]*rate_eval.p_p_to_d_weak_electron_capture"

    def test_function_string(self, rate1, rate2, rate3, rate4):

        ostr1 = """
@numba.njit()
def p_C13_to_N14(rate_eval, tf):
    # C13 + p --> N14
    rate = 0.0

    # nacrr
    rate += np.exp(  15.1825 + -13.5543*tf.T9i
                  + -1.5*tf.lnT9)
    # nacrn
    rate += np.exp(  18.5155 + -13.72*tf.T913i + -0.450018*tf.T913
                  + 3.70823*tf.T9 + -1.70545*tf.T953 + -0.666667*tf.lnT9)
    # nacrr
    rate += np.exp(  13.9637 + -5.78147*tf.T9i + -0.196703*tf.T913
                  + 0.142126*tf.T9 + -0.0238912*tf.T953 + -1.5*tf.lnT9)

    rate_eval.p_C13_to_N14 = rate
"""

        ostr2 = """
@numba.njit()
def p_p_He4_to_He3_He3(rate_eval, tf):
    # p + p + He4 --> He3 + He3
    rate = 0.0

    # nacrn
    rate += np.exp(  2.98257 + -149.222*tf.T9i + -12.277*tf.T913i + -0.103699*tf.T913
                  + -0.0649967*tf.T9 + 0.0168191*tf.T953 + -2.16667*tf.lnT9)

    rate_eval.p_p_He4_to_He3_He3 = rate
"""

        ostr3 = """
@numba.njit()
def n_p_He4_He4_to_He3_Li7(rate_eval, tf):
    # n + p + He4 + He4 --> He3 + Li7
    rate = 0.0

    # mafon
    rate += np.exp(  -14.8862 + -111.725*tf.T9i + -17.989*tf.T913i + -1.57523e-09*tf.T913
                  + 1.45934e-10*tf.T9 + -1.15341e-11*tf.T953 + -3.66667*tf.lnT9)

    rate_eval.n_p_He4_He4_to_He3_Li7 = rate
"""

        ostr4 = """
@numba.njit()
def p_p_to_d_weak_electron_capture(rate_eval, tf):
    # p + p --> d
    rate = 0.0

    # bet+w
    rate += np.exp(  -34.7863 + -3.51193*tf.T913i + 3.10086*tf.T913
                  + -0.198314*tf.T9 + 0.0126251*tf.T953 + -1.02517*tf.lnT9)
    #   ecw
    rate += np.exp(  -43.6499 + -0.00246064*tf.T9i + -2.7507*tf.T913i + -0.424877*tf.T913
                  + 0.015987*tf.T9 + -0.000690875*tf.T953 + -0.207625*tf.lnT9)

    rate_eval.p_p_to_d_weak_electron_capture = rate
"""

        assert rate1.function_string_py().strip() == ostr1.strip()
        assert rate2.function_string_py().strip() == ostr2.strip()
        assert rate3.function_string_py().strip() == ostr3.strip()
        assert rate4.function_string_py().strip() == ostr4.strip()
