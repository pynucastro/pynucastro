# unit tests for rates
import pytest
import warnings
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

    @pytest.fixture(scope='class')
    def pynet2(self, reaclib_library):

        n14agf18 = reaclib_library.get_rate_by_name("n14(a,g)f18")
        n14_new = pyna.ModifiedRate(n14agf18, new_products=["ne20"],
                                    stoichiometry={pyna.Nucleus("he4"): 1.5})
        new_n14_reverse = pyna.DerivedRate(rate=n14_new, compute_Q=True, use_pf=True)

        my_net = pyna.Library(rates=[n14_new, new_n14_reverse])
        pynet = pyna.PythonNetwork(libraries=my_net)
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

    def test_derived_modified_ydot_string(self, pynet2):
        ostr1 = \
"""dYdt[jhe4] = (
      ( + -1.5*rho*Y[jhe4]*Y[jn14]*rate_eval.He4_N14__Ne20__modified + 1.5*Y[jne20]*rate_eval.Ne20__He4_N14__derived )
   )

"""
        ostr2 = \
"""dYdt[jn14] = (
      ( -rho*Y[jhe4]*Y[jn14]*rate_eval.He4_N14__Ne20__modified +Y[jne20]*rate_eval.Ne20__He4_N14__derived )
   )

"""
        ostr3 = \
"""dYdt[jne20] = (
      ( +rho*Y[jhe4]*Y[jn14]*rate_eval.He4_N14__Ne20__modified -Y[jne20]*rate_eval.Ne20__He4_N14__derived )
   )

"""
        assert pynet2.full_ydot_string(pyna.Nucleus("he4"))  == ostr1
        assert pynet2.full_ydot_string(pyna.Nucleus("n14"))  == ostr2
        assert pynet2.full_ydot_string(pyna.Nucleus("ne20")) == ostr3

    def test_derived_modified_function_string(self, pynet2):

        ostr = \
"""@numba.njit()
def Ne20__He4_N14__derived(rate_eval, tf):
    # Ne20 --> 1.5 He4 + N14
    rate = 0.0

    # il10c
    rate += np.exp(  39.55827158733315 + -168.12237220574448*tf.T9i + -5.6227*tf.T913i)
    # il10c
    rate += np.exp(  25.85560958733315 + -162.31711220574448*tf.T9i)
    # il10c
    rate += np.exp(  47.19267158733315 + -157.1567722057445*tf.T9i + -36.2504*tf.T913i
                  + -5.0*tf.T953 + 0.833333*tf.lnT9)

    rate_eval.Ne20__He4_N14__derived = rate


    # setting He4 partition function to 1.0 by default, independent of T
    He4_pf = 1.0

    # interpolating Ne20 partition function
    Ne20_pf_exponent = np.interp(tf.T9, xp=Ne20_temp_array, fp=np.log10(Ne20_pf_array))
    Ne20_pf = 10.0**Ne20_pf_exponent

    # setting N14 partition function to 1.0 by default, independent of T
    N14_pf = 1.0

    z_r = He4_pf*N14_pf
    z_p = Ne20_pf
    rate_eval.Ne20__He4_N14__derived *= z_r/z_p
"""

        # To ignore partition function warning.
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=UserWarning)

            r = pynet2.get_rate_by_name("ne20(g,a)n14")
            assert r.function_string_py() == ostr
