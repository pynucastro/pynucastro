# unit tests for rates
import pytest


class TestPythonNetwork:
    @pytest.fixture(scope="class")
    @classmethod
    def rate1(cls, reaclib_library):
        return reaclib_library.get_rate_by_name("c13(p,g)n14")

    @pytest.fixture(scope="class")
    @classmethod
    def rate2(cls, reaclib_library):
        return reaclib_library.get_rate_by_name("he4(pp,he3)he3")

    @pytest.fixture(scope="class")
    @classmethod
    def rate3(cls, reaclib_library):
        return reaclib_library.get_rate_by_name("he4(npa,he3)li7")

    @pytest.fixture(scope="class")
    @classmethod
    def rate4(cls, reaclib_library):
        return reaclib_library.get_rate_by_name("p(p,)d")[1]

    def test_ydot_string(self, rate1, rate2, rate3, rate4):
        ydot1 = rate1.ydot_string_py()
        assert ydot1 == "rho*Y[jp]*Y[jc13]*rate_eval.C13_p_to_N14_reaclib"

        ydot2 = rate2.ydot_string_py()
        assert ydot2 == "5.00000000000000e-01*rho**2*Y[jp]**2*Y[jhe4]*rate_eval.He4_p_p_to_He3_He3_reaclib"

        ydot3 = rate3.ydot_string_py()
        assert ydot3 == "5.00000000000000e-01*rho**3*Y[jn]*Y[jp]*Y[jhe4]**2*rate_eval.He4_He4_p_n_to_He3_Li7_reaclib"

        ydot4 = rate4.ydot_string_py()
        assert ydot4 == "5.00000000000000e-01*rho**2*ye(Y)*Y[jp]**2*rate_eval.p_p_to_d_electron_capture_reaclib"

    def test_jacobian_string(self, rate1, rate2, rate3, rate4):
        jac1 = rate1.jacobian_string_py(rate1.reactants[0])
        assert jac1 == "rho*Y[jc13]*rate_eval.C13_p_to_N14_reaclib"

        jac2 = rate2.jacobian_string_py(rate2.reactants[0])
        assert jac2 == "5.00000000000000e-01*rho**2*2*Y[jp]*Y[jhe4]*rate_eval.He4_p_p_to_He3_He3_reaclib"

        jac3 = rate3.jacobian_string_py(rate3.reactants[0])
        assert jac3 == "5.00000000000000e-01*rho**3*Y[jp]*Y[jhe4]**2*rate_eval.He4_He4_p_n_to_He3_Li7_reaclib"

        jac4 = rate4.jacobian_string_py(rate4.reactants[0])
        assert jac4 == "5.00000000000000e-01*rho**2*ye(Y)*2*Y[jp]*rate_eval.p_p_to_d_electron_capture_reaclib"

    def test_function_string(self, rate1, rate2, rate3, rate4):

        ostr1 = """
@numba.njit()
def C13_p_to_N14_reaclib(rate_eval, tf, log_scor=0.0):
    # C13 + p --> N14
    rate = 0.0

    # nacrn
    ln_set_rate =  18.5155 + -13.72*tf.T913i + -0.450018*tf.T913 \\
                         + 3.70823*tf.T9 + -1.70545*tf.T953 + -0.666667*tf.lnT9

    ln_set_rate += log_scor
    set_rate = np.exp(ln_set_rate)
    rate += set_rate

    # nacrr
    ln_set_rate =  13.9637 + -5.78147*tf.T9i + -0.196703*tf.T913 \\
                         + 0.142126*tf.T9 + -0.0238912*tf.T953 + -1.5*tf.lnT9

    ln_set_rate += log_scor
    set_rate = np.exp(ln_set_rate)
    rate += set_rate

    # nacrr
    ln_set_rate =  15.1825 + -13.5543*tf.T9i \\
                         + -1.5*tf.lnT9

    ln_set_rate += log_scor
    set_rate = np.exp(ln_set_rate)
    rate += set_rate

    rate_eval.C13_p_to_N14_reaclib = rate
"""

        ostr2 = """
@numba.njit()
def He4_p_p_to_He3_He3_reaclib(rate_eval, tf, log_scor=0.0):
    # He4 + p + p --> He3 + He3
    rate = 0.0

    # nacrn
    ln_set_rate =  2.98257 + -149.222*tf.T9i + -12.277*tf.T913i + -0.103699*tf.T913 \\
                         + -0.0649967*tf.T9 + 0.0168191*tf.T953 + -2.16667*tf.lnT9

    ln_set_rate += log_scor
    set_rate = np.exp(ln_set_rate)
    rate += set_rate

    rate_eval.He4_p_p_to_He3_He3_reaclib = rate
"""

        ostr3 = """
@numba.njit()
def He4_He4_p_n_to_He3_Li7_reaclib(rate_eval, tf, log_scor=0.0):
    # He4 + He4 + p + n --> He3 + Li7
    rate = 0.0

    # mafon
    ln_set_rate =  -14.8862 + -111.725*tf.T9i + -17.989*tf.T913i + -1.57523e-09*tf.T913 \\
                         + 1.45934e-10*tf.T9 + -1.15341e-11*tf.T953 + -3.66667*tf.lnT9

    ln_set_rate += log_scor
    set_rate = np.exp(ln_set_rate)
    rate += set_rate

    rate_eval.He4_He4_p_n_to_He3_Li7_reaclib = rate
"""

        ostr4 = """
@numba.njit()
def p_p_to_d_electron_capture_reaclib(rate_eval, tf, log_scor=0.0):
    # p + p --> d
    rate = 0.0

    #   ecw
    ln_set_rate =  -43.6499 + -0.00246064*tf.T9i + -2.7507*tf.T913i + -0.424877*tf.T913 \\
                         + 0.015987*tf.T9 + -0.000690875*tf.T953 + -0.207625*tf.lnT9

    ln_set_rate += log_scor
    set_rate = np.exp(ln_set_rate)
    rate += set_rate

    rate_eval.p_p_to_d_electron_capture_reaclib = rate
"""

        assert rate1.function_string_py().strip() == ostr1.strip()
        assert rate2.function_string_py().strip() == ostr2.strip()
        assert rate3.function_string_py().strip() == ostr3.strip()
        assert rate4.function_string_py().strip() == ostr4.strip()
