# unit tests for rates
import pynucastro as pyna


class TestPythonDerivedNetwork:
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

        all_reactions = pyna.ReacLibLibrary()
        fwd_reactions = all_reactions.derived_forward()

        intermediate_nuclei = ['cr48', 'mn51', 'fe52', 'ni56', 'p', 'n', 'he4', 'co55']
        beta_decay_modes = ['fe54', 'fe55', 'fe54', 'fe56', 'cr56', 'mn56']

        all_nuclei = intermediate_nuclei + beta_decay_modes

        fwd_rates_lib = fwd_reactions.linking_nuclei(nuclist=all_nuclei, with_reverse=False)

        derived = []
        for r in fwd_rates_lib.get_rates():
            d = pyna.rates.DerivedRate(rate=r, use_A_nuc=False, use_pf=True)
            derived.append(d)

        der_rates_lib = pyna.Library(rates=derived)

        full_library = fwd_rates_lib + der_rates_lib

        pynet = pyna.PythonNetwork(libraries=full_library)
        self.pynet = pynet

    def teardown_method(self):
        """ this is run after each test """
        self.pynet = None

    def test_num_rates(self):
        assert len(self.pynet.rates) == 28

    def test_num_reaclib_rates(self):
        assert len(self.pynet.reaclib_rates) == 14

    def test_num_derived_rates(self):
        assert len(self.pynet.derived_rates) == 14

    def test_full_ydot_string(self):

        ostr = \
"""dYdt[jcr48] = (
   -rho*Y[jhe4]*Y[jcr48]*lambda_he4_cr48__fe52
   -rho*Y[jhe4]*Y[jcr48]*lambda_he4_cr48__p_mn51
   +Y[jfe52]*lambda_fe52__he4_cr48__derived
   +rho*Y[jp]*Y[jmn51]*lambda_p_mn51__he4_cr48__derived
   )

"""

        assert self.pynet.full_ydot_string(pyna.Nucleus("cr48")) == ostr

    def test_approx_function_string(self):

        ostr = \
"""@numba.njit()
def fe52__he4_cr48__derived(tf):
    # fe52 --> he4 + cr48
    rate = 0.0

    # ths8r
    rate += np.exp(  90.14738712482466 + -92.10912191363732*tf.T9i + -86.7459*tf.T913i + -9.79373*tf.T913
                  + -0.772169*tf.T9 + 0.155883*tf.T953 + 0.833333*tf.lnT9)

    he4_temp_array = np.array([0.01, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0])
    he4_pf_array = np.array([1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0])
    he4_pf_exponent = np.interp(tf.T9, xp=he4_temp_array, fp=np.log10(he4_pf_array))
    he4_pf = 10.0**he4_pf_exponent

    cr48_temp_array = np.array([0.01, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0, 22.0, 24.0, 26.0, 28.0, 30.0, 35.0, 40.0, 45.0, 50.0, 55.0, 60.0, 65.0, 70.0, 75.0, 80.0, 85.0, 90.0, 95.0, 100.0, 105.0, 110.0, 115.0, 120.0, 125.0, 130.0, 135.0, 140.0, 145.0, 150.0, 155.0, 160.0, 165.0, 170.0, 175.0, 180.0, 190.0, 200.0, 210.0, 220.0, 230.0, 240.0, 250.0, 275.0])
    cr48_pf_array = np.array([1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.000002, 1.000019, 1.000091, 1.000307, 1.000809, 1.01485, 1.063783, 1.15386, 1.279281, 1.432145, 1.60626, 1.797679, 2.004396, 2.46335, 2.995063, 3.632788, 4.440947, 5.538404, 9.62, 20.4, 50.9, 138.0, 378.0, 1030.0, 2740.0, 7160.0, 18400.0, 46400.0, 445000.0, 4040000.0, 35500000.0, 305000000.0, 2600000000.0, 21900000000.0, 184000000000.0, 1530000000000.0, 12700000000000.0, 105000000000000.0, 861000000000000.0, 7020000000000000.0, 5.69e+16, 4.59e+17, 3.68e+18, 2.94e+19, 2.34e+20, 1.85e+21, 1.45e+22, 1.14e+23, 8.9e+23, 6.93e+24, 5.38e+25, 4.16e+26, 3.22e+27, 2.48e+28, 1.91e+29, 1.46e+30, 1.12e+31, 8.6e+31, 5.03e+33, 2.94e+35, 1.71e+37, 9.98e+38, 5.82e+40, 3.41e+42, 2e+44, 5.34e+48])
    cr48_pf_exponent = np.interp(tf.T9, xp=cr48_temp_array, fp=np.log10(cr48_pf_array))
    cr48_pf = 10.0**cr48_pf_exponent

    fe52_temp_array = np.array([0.01, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0, 22.0, 24.0, 26.0, 28.0, 30.0, 35.0, 40.0, 45.0, 50.0, 55.0, 60.0, 65.0, 70.0, 75.0, 80.0, 85.0, 90.0, 95.0, 100.0, 105.0, 110.0, 115.0, 120.0, 125.0, 130.0, 135.0, 140.0, 145.0, 150.0, 155.0, 160.0, 165.0, 170.0, 175.0, 180.0, 190.0, 200.0, 210.0, 220.0, 230.0, 240.0, 250.0, 275.0])
    fe52_pf_array = np.array([1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.000004, 1.000022, 1.000087, 1.000261, 1.006988, 1.036149, 1.097022, 1.187924, 1.302823, 1.435985, 1.583461, 1.743116, 2.097506, 2.50807, 3.004926, 3.650892, 4.565183, 8.23, 18.6, 49.1, 139.0, 399.0, 1140.0, 3190.0, 8840.0, 24200.0, 65800.0, 774000.0, 8800000.0, 97700000.0, 1070000000.0, 11500000000.0, 122000000000.0, 1290000000000.0, 13400000000000.0, 138000000000000.0, 1410000000000000.0, 1.42e+16, 1.42e+17, 1.41e+18, 1.39e+19, 1.36e+20, 1.32e+21, 1.27e+22, 1.22e+23, 1.16e+24, 1.1e+25, 1.04e+26, 9.76e+26, 9.14e+27, 8.53e+28, 7.93e+29, 7.36e+30, 6.82e+31, 6.3e+32, 5.81e+33, 5.35e+34, 4.52e+36, 3.81e+38, 3.2e+40, 2.69e+42, 2.26e+44, 1.91e+46, 1.61e+48, 1.07e+53])
    fe52_pf_exponent = np.interp(tf.T9, xp=fe52_temp_array, fp=np.log10(fe52_pf_array))
    fe52_pf = 10.0**fe52_pf_exponent


    z_r = he4_pf*cr48_pf
    z_p = fe52_pf
    rate *= z_r/z_p

    return rate

"""

        r = self.pynet.get_rate("fe52__he4_cr48__derived")
        assert r.function_string_py() == ostr
