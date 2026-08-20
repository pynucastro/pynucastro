# tests of network containing the rate he4(npa,p)be9 and its inverse
# since protons are on both sides, the contribution to dY(p)/dt should
# be zero.  We assert that here

import pytest

import pynucastro as pyna


class TestRateNetZero:

    @pytest.fixture(scope="class")
    @classmethod
    def lib(cls, reaclib_library):
        rr = reaclib_library.get_rate_by_name("he4(npa,p)be9")
        rv = reaclib_library.get_rate_by_name("be9(p,npa)he4")

        return pyna.Library(rates=[rr, rv])

    def test_python_net(self, lib):
        net = pyna.PythonNetwork(libraries=lib)

        # protons should cancel out in their contribution
        Ydot_p_str = net.full_ydot_string(pyna.Nucleus("p"))

        assert Ydot_p_str == "dYdt[jp] = (\n      0.0\n   )\n\n"

    def test_cxx_net(self, lib):
        net = pyna.SimpleCxxNetwork(libraries=lib)

        net.compose_ydot()
        assert net.ydot_out_result[pyna.Nucleus("p")] == [(None, None)]

