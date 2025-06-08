# unit tests for rates
import pytest

import pynucastro as pyna


class TestPythonDerivedNetwork:
    @pytest.fixture(scope='class')
    def pynet(self, reaclib_library):

        fwd_reactions = reaclib_library.forward_for_detailed_balance()

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
      ( -rho*Y[jhe4]*Y[jcr48]*rate_eval.He4_Cr48__Fe52 +Y[jfe52]*rate_eval.Fe52__He4_Cr48__derived ) +
      ( -rho*Y[jhe4]*Y[jcr48]*rate_eval.He4_Cr48__p_Mn51 +rho*Y[jp]*Y[jmn51]*rate_eval.p_Mn51__He4_Cr48__derived )
   )

"""

        assert pynet.full_ydot_string(pyna.Nucleus("cr48")) == ostr

    def test_derived_function_string(self, pynet):

        ostr = \
"""@numba.njit()
def Fe52__p_Mn51__derived(rate_eval, tf):
    # Fe52 --> p + Mn51
    rate = 0.0

    # ths8r
    rate += np.exp(  61.74743132228039 + -85.63264034844842*tf.T9i + -36.1825*tf.T913i + 0.873042*tf.T913
                  + -2.89731*tf.T9 + 0.364394*tf.T953 + 0.833333*tf.lnT9)

    rate_eval.Fe52__p_Mn51__derived = rate


    # interpolating Mn51 partition function
    Mn51_pf_exponent = np.interp(tf.T9, xp=Mn51_temp_array, fp=np.log10(Mn51_pf_array))
    Mn51_pf = 10.0**Mn51_pf_exponent

    # setting p partition function to 1.0 by default, independent of T
    p_pf = 1.0

    # interpolating Fe52 partition function
    Fe52_pf_exponent = np.interp(tf.T9, xp=Fe52_temp_array, fp=np.log10(Fe52_pf_array))
    Fe52_pf = 10.0**Fe52_pf_exponent

    z_r = p_pf*Mn51_pf
    z_p = Fe52_pf
    rate_eval.Fe52__p_Mn51__derived *= z_r/z_p
"""

        r = pynet.get_rate("fe52__p_mn51__derived")
        assert r.function_string_py() == ostr
