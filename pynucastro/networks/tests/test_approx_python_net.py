# unit tests for rates
import pytest

import pynucastro as pyna


class TestPythonNetwork:
    @pytest.fixture(scope="class")
    def pynet(self, reaclib_library):
        mynet = reaclib_library.linking_nuclei(["p", "he4", "mg24",
                                                "al27", "si28", "p31", "s32"])
        pynet = pyna.PythonNetwork(libraries=[mynet])
        pynet.make_ap_pg_approx()
        pynet.remove_nuclei(["al27", "p31"])
        return pynet

    def test_num_rates(self, pynet):
        assert len(pynet.rates) == 4

    def test_num_reaclib_rates(self, pynet):
        assert len(pynet.reaclib_rates) == 12

    def test_num_approx_rates(self, pynet):
        assert len(pynet.approx_rates) == 4

    def test_full_ydot_string(self, pynet):
        ostr = \
"""dYdt[jhe4] = (
   -rho*Y[jhe4]*Y[jmg24]*rate_eval.mg24_he4__si28__approx
   -rho*Y[jhe4]*Y[jsi28]*rate_eval.si28_he4__s32__approx
   +Y[jsi28]*rate_eval.si28__mg24_he4__approx
   +Y[js32]*rate_eval.s32__si28_he4__approx
   )

"""

        assert pynet.full_ydot_string(pyna.Nucleus("he4")) == ostr

    def test_approx_function_string(self, pynet):

        ostr = \
"""@numba.njit()
def mg24_he4__si28__approx(rate_eval, tf):
    r_ag = rate_eval.he4_mg24__si28
    r_ap = rate_eval.he4_mg24__p_al27
    r_pg = rate_eval.p_al27__si28
    r_pa = rate_eval.p_al27__he4_mg24
    rate = r_ag + r_ap * r_pg / (r_pg + r_pa)
    rate_eval.mg24_he4__si28__approx = rate

"""
        r = pynet.get_rate("mg24_he4__si28__approx")
        assert r.function_string_py() == ostr

    def test_function_string(self, pynet):

        ostr = \
"""@numba.njit()
def he4_mg24__si28(rate_eval, tf):
    # mg24 + he4 --> si28
    rate = 0.0

    # st08r
    rate += np.exp(  -50.5494 + -12.8332*tf.T9i + 21.3721*tf.T913i + 37.7649*tf.T913
                  + -4.10635*tf.T9 + 0.249618*tf.T953 + -1.5*tf.lnT9)
    # st08r
    rate += np.exp(  8.03977 + -15.629*tf.T9i
                  + -1.5*tf.lnT9)

    rate_eval.he4_mg24__si28 = rate

"""

        r = pynet.get_rate("mg24_he4__si28__approx")
        assert r.get_child_rates()[0].function_string_py().strip() == ostr.strip()
