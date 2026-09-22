# unit tests for rates
import importlib
import sys
from pathlib import Path

import numpy as np
import pytest
from pytest import approx
from scipy.integrate import solve_ivp

import pynucastro as pyna
from pynucastro.screening import chugunov_2007


@pytest.mark.parametrize("direction", ["forward", "reverse", "both"])
def test_double_neutron_capture_jacobian(reaclib_library, tmp_path, direction):
    """Check dYdot/dY, including implicit neutron dependence, at fixed rho/T."""
    library = reaclib_library.linking_nuclei(["n", "fe52", "fe53", "fe54"])
    net = pyna.PythonNetwork(libraries=[library])
    net.make_nn_g_approx(intermediate_nuclei=["fe53"])
    assert len(net.approx_rates) == 2
    if direction != "both":
        rates = [r for r in net.approx_rates if r.is_reverse == (direction == "reverse")]
        net = pyna.PythonNetwork(rates=rates)

    path = tmp_path / "nn_capture.py"
    net.write_network(path)
    spec = importlib.util.spec_from_file_location("nn_capture", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)

    rho = 1.e7
    temperature = 3.e9
    abundances = np.array([0.2 if n.raw == "n" else 0.4 / n.A for n in net.unique_nuclei])
    # Evaluate the generated analytic Jacobian, including the product-rule
    # terms from the approximate rates' implicit neutron dependence.
    jac = module.jacobian(0.0, abundances, rho, temperature)
    for j, abundance in enumerate(abundances):
        # Perturb only Y_j at fixed density and temperature. The centered
        # finite difference of the RHS approximates column j of dYdot/dY.
        step = 1.e-5 * abundance
        plus = abundances.copy()
        minus = abundances.copy()
        plus[j] += step
        minus[j] -= step
        # Each RHS evaluation recomputes the composition-dependent rates,
        # so this independently checks their implicit abundance derivatives.
        numerical = (module.rhs(0.0, plus, rho, temperature) -
                     module.rhs(0.0, minus, rho, temperature)) / (2 * step)
        # Compare the finite-difference column with the analytic column.
        np.testing.assert_allclose(jac[:, j], numerical, rtol=1.e-7, atol=1.e-8)


class TestPythonNetwork:
    @pytest.fixture(scope="class")
    @classmethod
    def pynet(cls, reaclib_library):
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
      ( -rho*Y[jhe4]*Y[jmg24]*rate_eval.Mg24_He4_to_Si28_approx +Y[jsi28]*rate_eval.Si28_to_Mg24_He4_approx ) +
      ( -rho*Y[jhe4]*Y[jsi28]*rate_eval.Si28_He4_to_S32_approx +Y[js32]*rate_eval.S32_to_Si28_He4_approx )
   )

"""

        assert pynet.full_ydot_string(pyna.Nucleus("he4")) == ostr

    def test_approx_function_string(self, pynet):

        ostr = \
"""@numba.njit()
def Mg24_He4_to_Si28_approx(rate_eval):
    r_pg = rate_eval.Al27_p_to_Si28_reaclib
    r_pa = rate_eval.Al27_p_to_He4_Mg24_reaclib
    r_pY = 0.0
    r_ag = rate_eval.Mg24_He4_to_Si28_reaclib
    r_ap = rate_eval.Mg24_He4_to_p_Al27_reaclib
    rate = r_ag + r_ap * r_pg / (r_pg + r_pa + r_pY)
    rate_eval.Mg24_He4_to_Si28_approx = rate

"""
        r = pynet.get_rate("mg24_he4_to_si28_approx")
        assert r.function_string_py() == ostr

    def test_function_string(self, pynet):

        ostr = \
"""@numba.njit()
def Mg24_He4_to_Si28_reaclib(rate_eval, tf, log_scor=0.0):
    # Mg24 + He4 --> Si28
    rate = 0.0

    # st08r
    ln_set_rate =  8.03977 + -15.629*tf.T9i \\
                         + -1.5*tf.lnT9

    ln_set_rate += log_scor
    set_rate = np.exp(ln_set_rate)
    rate += set_rate

    # st08r
    ln_set_rate =  -50.5494 + -12.8332*tf.T9i + 21.3721*tf.T913i + 37.7649*tf.T913 \\
                         + -4.10635*tf.T9 + 0.249618*tf.T953 + -1.5*tf.lnT9

    ln_set_rate += log_scor
    set_rate = np.exp(ln_set_rate)
    rate += set_rate

    rate_eval.Mg24_He4_to_Si28_reaclib = rate

"""

        r = pynet.get_rate("mg24_he4_to_si28_approx")
        print(r)
        assert r.get_child_rates()[1].function_string_py().strip() == ostr.strip()

    def test_integrating(self, pynet):
        pynet.write_network("app.py")
        app = importlib.import_module("app")

        rho = 1.e7
        T = 3e9

        X0 = np.zeros(app.nnuc)
        X0[app.jhe4] = 0.5
        X0[app.jmg24] = 0.5

        Y0 = X0 / app.A

        tmax = 1.e-3
        sol = solve_ivp(app.rhs, [0, tmax], Y0, method="BDF",
                        jac=app.jacobian,
                        dense_output=True, args=(rho, T, chugunov_2007), rtol=1.e-6, atol=1.e-10)

        # these are the final molar fractions
        answer = [8.33333490e-02, 9.24569852e-20, 1.56798113e-08, 2.08333177e-02]

        for i in range(app.nnuc):
            assert answer[i] == approx(sol.y[i, -1])

        # clean up generated files if the test passed
        Path("app.py").unlink()
        # remove imported module from cache
        del app
        del sys.modules["app"]

    def test_to_composition(self, pynet):
        pynet.write_network("app2.py")
        app2 = importlib.import_module("app2")

        comp_orig = pyna.Composition(pynet.unique_nuclei)
        comp_orig.set_solar_like()

        Y = np.zeros(app2.nnuc)
        for nuc, molar_fraction in comp_orig.get_molar().items():
            Y[app2.names.index(nuc.caps_name)] = molar_fraction
        comp_new = app2.to_composition(Y)

        assert comp_new.X == comp_orig.X

        # clean up generated files if the test passed
        Path("app2.py").unlink()
        # remove imported module from cache
        del app2
        del sys.modules["app2"]
