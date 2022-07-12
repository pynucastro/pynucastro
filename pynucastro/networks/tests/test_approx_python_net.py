# unit tests for rates
import pynucastro as pyna


class TestPythonNetwork:
    @classmethod
    def setup_class(cls):
        """ this is run once for each class before any tests """
        pass

    @classmethod
    def teardown_class(cls):
        """ this is run once for each class after all tests """
        pass

    def setup_method(self):
        """ this is run before each test """

        rl = pyna.ReacLibLibrary()
        mynet = rl.linking_nuclei(["p", "he4", "mg24",
                                   "al27", "si28", "p31", "s32"])
        pynet = pyna.PythonNetwork(libraries=[mynet])
        pynet.make_ap_pg_approx()
        pynet.remove_nuclei(["al27", "p31"])
        self.pynet = pynet

    def teardown_method(self):
        """ this is run after each test """
        self.pynet = None

    def test_num_rates(self):
        assert len(self.pynet.rates) == 4

    def test_full_ydot_string(self):
        ostr = \
"""dYdt[jhe4] = (
   -rho*Y[jhe4]*Y[jmg24]*lambda_mg24_he4__si28__approx
   -rho*Y[jhe4]*Y[jsi28]*lambda_si28_he4__s32__approx
   +Y[jsi28]*lambda_si28__mg24_he4__approx
   +Y[js32]*lambda_s32__si28_he4__approx
   )

"""

        assert self.pynet.full_ydot_string(pyna.Nucleus("he4")) == ostr

    def test_approx_function_string(self):

        ostr = \
"""@numba.njit()
def mg24_he4__si28__approx(tf):
    r_ag = he4_mg24__si28(tf)
    r_ap = he4_mg24__p_al27(tf)
    r_pg = p_al27__si28(tf)
    r_pa = p_al27__he4_mg24(tf)
    rate = r_ag + r_ap * r_pg / (r_pg + r_pa)
    return rate

"""
        r = self.pynet.get_rate("mg24_he4__si28__approx")
        assert r.function_string_py() == ostr

    def test_function_string(self):

        ostr = \
"""@numba.njit()
def he4_mg24__si28(tf):
    # mg24 + he4 --> si28
    rate = 0.0

    # st08r
    rate += np.exp(  -50.5494 + -12.8332*tf.T9i + 21.3721*tf.T913i + 37.7649*tf.T913
                  + -4.10635*tf.T9 + 0.249618*tf.T953 + -1.5*tf.lnT9)
    # st08r
    rate += np.exp(  8.03977 + -15.629*tf.T9i
                  + -1.5*tf.lnT9)

    return rate

"""

        r = self.pynet.get_rate("mg24_he4__si28__approx")
        assert r.get_child_rates()[0].function_string_py().strip() == ostr.strip()
