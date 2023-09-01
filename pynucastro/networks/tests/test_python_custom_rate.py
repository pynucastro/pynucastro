# unit tests for rates
import pytest
from pytest import approx

import pynucastro as pyna


class MyRate(pyna.Rate):
    def __init__(self, reactants=None, products=None,
                 r0=1.0, T0=1.0, nu=0):

        # we set the chapter to custom so the network knows how to deal with it
        self.chapter = "custom"

        # call the Rate init to do the remaining initialization
        super().__init__(reactants=reactants, products=products)

        self.r0 = r0
        self.T0 = T0
        self.nu = nu

    def function_string_py(self):
        """return a string containing a python function that computes
        the rate"""
        fstring = ""
        fstring += "@numba.njit()\n"
        fstring += f"def {self.fname}(rate_eval, tf):\n"
        fstring += f"    rate_eval.{self.fname} = {self.r0} * (tf.T9 * 1.e9 / {self.T0} )**({self.nu})\n\n"
        return fstring

    def eval(self, T, rhoY=None):
        return self.r0 * (T / self.T0)**self.nu


class TestPythonCustomNetwork:
    @pytest.fixture(scope="class")
    def pynet(self, reaclib_library):

        # create a CNO like network but with a custom rate
        # for N14(p,g)O15

        rate_names = ["c12(p,g)n13",
                      "c13(p,g)n14",
                      "n13(,)c13",
                      "n13(p,g)o14",
                      "n15(p,a)c12",
                      "o14(,)n14",
                      "o15(,)n15"]

        rates = reaclib_library.get_rate_by_name(rate_names)

        r = reaclib_library.get_rate_by_name("n14(p,g)o15")
        T0 = 3.e7
        nu = r.get_rate_exponent(T0)
        r0 = r.eval(T0)

        r_custom = MyRate(reactants=[pyna.Nucleus("n14"), pyna.Nucleus("p")],
                          products=[pyna.Nucleus("o15")],
                          r0=r0, T0=T0, nu=nu)

        return pyna.PythonNetwork(rates=rates+[r_custom])

    def test_eval(self, pynet):
        r_custom = pynet.get_rate_by_name("n14(p,g)o15")
        assert r_custom.eval(1.e8) == approx(2.0377211133509627e-5)

    def test_ydot_string(self, pynet):

        dpdt = \
"""dYdt[jp] = (
   -rho*Y[jp]*Y[jc12]*rate_eval.p_c12__n13
   -rho*Y[jp]*Y[jc13]*rate_eval.p_c13__n14
   -rho*Y[jp]*Y[jn13]*rate_eval.p_n13__o14
   -rho*Y[jp]*Y[jn15]*rate_eval.p_n15__he4_c12
   -rho*Y[jp]*Y[jn14]*rate_eval.n14_p__o15__generic
   )

"""

        assert pynet.full_ydot_string(pyna.Nucleus("p")) == dpdt

    def test_py_function(self, pynet):
        r_custom = pynet.get_rate_by_name("n14(p,g)o15")

        func = \
"""@numba.njit()
def n14_p__o15__generic(rate_eval, tf):
    rate_eval.n14_p__o15__generic = 1.416655077954945e-13 * (tf.T9 * 1.e9 / 30000000.0 )**(15.601859314950396)

"""

        assert r_custom.function_string_py() == func
