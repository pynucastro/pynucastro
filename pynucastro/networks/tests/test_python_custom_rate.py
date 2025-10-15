# unit tests for rates
import importlib
import sys
from pathlib import Path

import numpy as np
import pytest
from pytest import approx

import pynucastro as pyna
from pynucastro.screening import chugunov_2007

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

    def eval(self, T, *, rho=None, comp=None,
             screen_func=None, symmetric_screening=False):
        """Evaluate the rate along with screening correction."""

        r = self.r0 * (T / self.T0)**self.nu

        scor = 1.0
        if screen_func is not None:
            if rho is None or comp is None:
                raise ValueError("rho (density) and comp (Composition) needs to be defined when applying electron screening.")
            scor = self.evaluate_screening(rho, T, comp, screen_func,
                                           symmetric_screening=symmetric_screening)
        r *= scor
        return r


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

    def test_eval_screening(self, pynet):
        T = 1.e8
        rho = 1.e7
        comp = pyna.Composition(pynet.unique_nuclei)
        comp.set_equal()

        r_custom = pynet.get_rate_by_name("n14(p,g)o15")
        r = r_custom.eval(T, rho=rho, comp=comp,
                          screen_func=chugunov_2007)

        assert r  == approx(2.0377211133509627e-5)

    def test_ydot_string(self, pynet):

        dpdt = \
"""dYdt[jp] = (
      -rho*Y[jp]*Y[jc12]*rate_eval.p_C12__N13  +
      -rho*Y[jp]*Y[jc13]*rate_eval.p_C13__N14  +
      -rho*Y[jp]*Y[jn13]*rate_eval.p_N13__O14  +
      -rho*Y[jp]*Y[jn15]*rate_eval.p_N15__He4_C12  +
      -rho*Y[jp]*Y[jn14]*rate_eval.N14_p__O15__generic
   )

"""

        assert pynet.full_ydot_string(pyna.Nucleus("p")) == dpdt

    def test_py_function(self, pynet):
        r_custom = pynet.get_rate_by_name("n14(p,g)o15")

        func = \
"""@numba.njit()
def N14_p__O15__generic(rate_eval, tf):
    rate_eval.N14_p__O15__generic = 1.416655077954945e-13 * (tf.T9 * 1.e9 / 30000000.0 )**(15.601859314950396)

"""

        assert r_custom.function_string_py() == func

    def test_evaluate_ydots(self, pynet):
        rho = 1.e4
        T = 4e7
        comp = pyna.Composition(pynet.unique_nuclei)
        comp.set_equal()

        ydots = pynet.evaluate_ydots(rho, T, comp)

        ydots_ref = {
            pyna.Nucleus("p"): -4.89131088e-06,
            pyna.Nucleus("he4"): 4.84655620e-06,
            pyna.Nucleus("c12"): 4.83556908e-06,
            pyna.Nucleus("c13"): 9.87368448e-06,
            pyna.Nucleus("n13"): -9.89635377e-06,
            pyna.Nucleus("n14"): 7.79536222e-05,
            pyna.Nucleus("n15"): 3.72390497e-05,
            pyna.Nucleus("o14"): -7.79200769e-05,
            pyna.Nucleus("o15"): -4.20854947e-05,
        }

        for n, value in ydots_ref.items():
            assert ydots[n] == approx(value)

    def test_import(self, pynet):
        pynet.write_network("custom.py")
        custom = importlib.import_module("custom")

        rho = 1.e4
        T = 4e7

        X0 = np.ones(custom.nnuc) / custom.nnuc
        Y0 = X0 / custom.A

        ydots = custom.rhs(0.0, Y0, rho, T)

        ydot_save = np.array([-4.89131088e-06, 4.84655620e-06,
                              4.83556908e-06, 9.87368448e-06,
                              -9.89635377e-06, 7.79536222e-05,
                              3.72390497e-05, -7.79200769e-05,
                              -4.20854947e-05])

        for n in range(custom.nnuc):
            assert ydots[n] == approx(ydot_save[n])

        Path("custom.py").unlink()
        del custom
        del sys.modules["custom"]
