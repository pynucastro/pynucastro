# unit tests for rates
import pynucastro as pyna
import numpy as np
import pytest


class TestPythonDerivedNetwork:
    @pytest.fixture(scope='class')
    def pynet(self, reaclib_library):

        fwd_reactions = reaclib_library.derived_forward()

        intermediate_nuclei = ['cr48', 'mn51', 'fe52', 'ni56', 'p', 'n', 'he4', 'co55']
        beta_decay_modes = ['fe54', 'fe55', 'fe54', 'fe56', 'cr56', 'mn56']

        all_nuclei = intermediate_nuclei + beta_decay_modes

        fwd_rates_lib = fwd_reactions.linking_nuclei(nuclist=all_nuclei, with_reverse=False)

        derived = []
        for r in fwd_rates_lib.get_rates():
            d = pyna.rates.DerivedRate(rate=r, compute_Q=False, use_pf=True)
            derived.append(d)

        der_rates_lib = pyna.Library(rates=derived)

        full_library = fwd_rates_lib + der_rates_lib

        pynet = pyna.PythonNetwork(libraries=full_library)
        return pynet

    def test_num_rates(self, pynet):
        assert len(pynet.rates) == 28

    def test_num_reaclib_rates(self, pynet):
        assert len(pynet.reaclib_rates) == 14

    def test_num_derived_rates(self, pynet):
        assert len(pynet.derived_rates) == 14

    def test_full_ydot_string(self, pynet):

        ostr = \
"""dYdt[jcr48] = (
   -rho*Y[jhe4]*Y[jcr48]*rate_eval.he4_cr48__fe52
   -rho*Y[jhe4]*Y[jcr48]*rate_eval.he4_cr48__p_mn51
   +Y[jfe52]*rate_eval.fe52__he4_cr48__derived
   +rho*Y[jp]*Y[jmn51]*rate_eval.p_mn51__he4_cr48__derived
   )

"""

        assert pynet.full_ydot_string(pyna.Nucleus("cr48")) == ostr

    def test_approx_function_string(self, pynet):

        dispute_values = np.array([61.72798916565748, 61.72798916565747])
        ostr = [""" """, """ """]
        for i in range(2):
            ostr[i] = \
"""@numba.njit()
def fe52__p_mn51__derived(rate_eval, tf):
    # fe52 --> p + mn51
    rate = 0.0

    # ths8r
    rate += np.exp(  {:.14f} + -85.6326403498911*tf.T9i + -36.1825*tf.T913i + 0.873042*tf.T913
                  + -2.89731*tf.T9 + 0.364394*tf.T953 + 0.833333*tf.lnT9)

    #Setting p partition function to 1.0 by default, independent of T
    p_pf = 1.0

    #Interpolating mn51 partition function
    mn51_temp_array = np.array([0.01, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0, 22.0, 24.0, 26.0, 28.0, 30.0, 35.0, 40.0, 45.0, 50.0, 55.0, 60.0, 65.0, 70.0, 75.0, 80.0, 85.0, 90.0, 95.0, 100.0, 105.0, 110.0, 115.0, 120.0, 125.0, 130.0, 135.0, 140.0, 145.0, 150.0, 155.0, 160.0, 165.0, 170.0, 175.0, 180.0, 190.0, 200.0, 210.0, 220.0, 230.0, 240.0, 250.0, 275.0])
    mn51_pf_array = np.array([1.0, 1.0, 1.000001, 1.000139, 1.001372, 1.005432, 1.013593, 1.026171, 1.042777, 1.062688, 1.085109, 1.213221, 1.339494, 1.454387, 1.561577, 1.666693, 1.774685, 1.889684, 2.015407, 2.314315, 2.707743, 3.248632, 4.018887, 5.147025, 9.37, 19.4, 43.7, 104.0, 253.0, 623.0, 1550.0, 3850.0, 9610.0, 24000.0, 236000.0, 2310000.0, 22600000.0, 220000000.0, 2130000000.0, 20500000000.0, 196000000000.0, 1860000000000.0, 17600000000000.0, 164000000000000.0, 1530000000000000.0, 1.41e+16, 1.29e+17, 1.18e+18, 1.07e+19, 9.67e+19, 8.68e+20, 7.76e+21, 6.91e+22, 6.13e+23, 5.42e+24, 4.78e+25, 4.2e+26, 3.69e+27, 3.23e+28, 2.83e+29, 2.47e+30, 2.16e+31, 1.88e+32, 1.64e+33, 1.24e+35, 9.43e+36, 7.17e+38, 5.46e+40, 4.17e+42, 3.2e+44, 2.47e+46, 1.32e+51])
    mn51_pf_exponent = np.interp(tf.T9, xp=mn51_temp_array, fp=np.log10(mn51_pf_array))
    mn51_pf = 10.0**mn51_pf_exponent

    #Interpolating fe52 partition function
    fe52_temp_array = np.array([0.01, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0, 22.0, 24.0, 26.0, 28.0, 30.0, 35.0, 40.0, 45.0, 50.0, 55.0, 60.0, 65.0, 70.0, 75.0, 80.0, 85.0, 90.0, 95.0, 100.0, 105.0, 110.0, 115.0, 120.0, 125.0, 130.0, 135.0, 140.0, 145.0, 150.0, 155.0, 160.0, 165.0, 170.0, 175.0, 180.0, 190.0, 200.0, 210.0, 220.0, 230.0, 240.0, 250.0, 275.0])
    fe52_pf_array = np.array([1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.000004, 1.000022, 1.000087, 1.000261, 1.006988, 1.036149, 1.097022, 1.187924, 1.302823, 1.435985, 1.583461, 1.743116, 2.097506, 2.50807, 3.004926, 3.650892, 4.565183, 8.23, 18.6, 49.1, 139.0, 399.0, 1140.0, 3190.0, 8840.0, 24200.0, 65800.0, 774000.0, 8800000.0, 97700000.0, 1070000000.0, 11500000000.0, 122000000000.0, 1290000000000.0, 13400000000000.0, 138000000000000.0, 1410000000000000.0, 1.42e+16, 1.42e+17, 1.41e+18, 1.39e+19, 1.36e+20, 1.32e+21, 1.27e+22, 1.22e+23, 1.16e+24, 1.1e+25, 1.04e+26, 9.76e+26, 9.14e+27, 8.53e+28, 7.93e+29, 7.36e+30, 6.82e+31, 6.3e+32, 5.81e+33, 5.35e+34, 4.52e+36, 3.81e+38, 3.2e+40, 2.69e+42, 2.26e+44, 1.91e+46, 1.61e+48, 1.07e+53])
    fe52_pf_exponent = np.interp(tf.T9, xp=fe52_temp_array, fp=np.log10(fe52_pf_array))
    fe52_pf = 10.0**fe52_pf_exponent

    z_r = p_pf*mn51_pf
    z_p = fe52_pf
    rate *= z_r/z_p

    rate_eval.fe52__p_mn51__derived = rate

""".format(dispute_values[i])

        r = pynet.get_rate("fe52__p_mn51__derived")
        assert r.function_string_py() == ostr[0] or r.function_string_py() == ostr[1]
