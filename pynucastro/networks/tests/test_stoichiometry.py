import io

import pytest

import pynucastro as pyna


class TestStoichiometry:
    """test the ydot generation when we explicitly change the
    stichiometry coefficients of a rate in a network"""

    @pytest.fixture(scope="class")
    def lib(self, reaclib_library):
        nuclei = ["he4", "c12", "o16"]
        lib = reaclib_library.linking_nuclei(nuclei, with_reverse=False)
        c12ag = lib.get_rate_by_name("c12(a,g)o16")
        c12ag.stoichiometry = {pyna.Nucleus("he4"): 10,
                               pyna.Nucleus("c12"): 20,
                               pyna.Nucleus("o16"): 40}
        return lib

    def test_python_ydot(self, lib):
        net = pyna.PythonNetwork(libraries=lib)

        dyhe4_dt = net.full_ydot_string(pyna.Nucleus("he4"))
        assert dyhe4_dt == """dYdt[jhe4] = (\n      + -10*rho*Y[jhe4]*Y[jc12]*rate_eval.He4_C12__O16  +\n      + -3*1.66666666666667e-01*rho**2*Y[jhe4]**3*rate_eval.He4_He4_He4__C12\n   )\n\n"""

        dyo16_dt = net.full_ydot_string(pyna.Nucleus("o16"))
        assert dyo16_dt == """dYdt[jo16] = (\n      + 40*rho*Y[jhe4]*Y[jc12]*rate_eval.He4_C12__O16\n   )\n\n"""

    def test_cxx_ydot(self, lib):
        net = pyna.SimpleCxxNetwork(libraries=lib)
        net.compose_ydot()

        output = io.StringIO()
        net._write_ydot_nuc(0, output, net.ydot_out_result[pyna.Nucleus("he4")])  # pylint: disable=protected-access
        result = output.getvalue()
        output.close()
        assert result == """-10.0*screened_rates(k_He4_C12_to_O16)*Y(C12)*Y(He4)*state.rho +\n-0.5*screened_rates(k_He4_He4_He4_to_C12)*std::pow(Y(He4), 3)*std::pow(state.rho, 2);\n\n"""

        output = io.StringIO()
        net._write_ydot_nuc(0, output, net.ydot_out_result[pyna.Nucleus("o16")])  # pylint: disable=protected-access
        result = output.getvalue()
        output.close()
        assert result == """40.0*screened_rates(k_He4_C12_to_O16)*Y(C12)*Y(He4)*state.rho;\n\n"""
