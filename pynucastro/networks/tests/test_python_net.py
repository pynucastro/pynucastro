# unit tests for rates
import pytest

from pynucastro import rates


class TestPythonNetwork:
    @pytest.fixture(scope="class")
    def rate(self):
        return rates.ReacLibRate("c13-pg-n14-nacr")

    def test_ydot_string(self, rate):
        ydot = rate.ydot_string_py()
        assert ydot in ("rho*Y[jc13]*Y[jp]*rate_eval.p_c13__n14",
                        "rho*Y[jp]*Y[jc13]*rate_eval.p_c13__n14")

    def test_jacobian_string(self, rate):
        jac = rate.jacobian_string_py(rate.reactants[0])
        assert jac == "rho*Y[jc13]*rate_eval.p_c13__n14"

    def test_function_string(self, rate):

        ostr = """
@numba.njit()
def p_c13__n14(rate_eval, tf):
    # c13 + p --> n14
    rate = 0.0

    # nacrn
    rate += np.exp(  18.5155 + -13.72*tf.T913i + -0.450018*tf.T913
                  + 3.70823*tf.T9 + -1.70545*tf.T953 + -0.666667*tf.lnT9)
    # nacrr
    rate += np.exp(  13.9637 + -5.78147*tf.T9i + -0.196703*tf.T913
                  + 0.142126*tf.T9 + -0.0238912*tf.T953 + -1.5*tf.lnT9)
    # nacrr
    rate += np.exp(  15.1825 + -13.5543*tf.T9i
                  + -1.5*tf.lnT9)

    rate_eval.p_c13__n14 = rate
"""

        assert rate.function_string_py().strip() == ostr.strip()
