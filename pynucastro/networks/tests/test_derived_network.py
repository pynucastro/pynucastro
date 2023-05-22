# unit tests for rates
import numpy as np
import pytest

import pynucastro as pyna


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

    rate_eval.fe52__p_mn51__derived = rate


    # interpolating mn51 partition function
    mn51_pf_exponent = np.interp(tf.T9, xp=mn51_temp_array, fp=np.log10(mn51_pf_array))
    mn51_pf = 10.0**mn51_pf_exponent

    # setting p partition function to 1.0 by default, independent of T
    p_pf = 1.0

    # interpolating fe52 partition function
    fe52_pf_exponent = np.interp(tf.T9, xp=fe52_temp_array, fp=np.log10(fe52_pf_array))
    fe52_pf = 10.0**fe52_pf_exponent

    z_r = p_pf*mn51_pf
    z_p = fe52_pf
    rate_eval.fe52__p_mn51__derived *= z_r/z_p
""".format(dispute_values[i])

        r = pynet.get_rate("fe52__p_mn51__derived")
        assert r.function_string_py() == ostr[0] or r.function_string_py() == ostr[1]
