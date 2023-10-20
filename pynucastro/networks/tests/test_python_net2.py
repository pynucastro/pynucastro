# unit tests for rates
import pytest

import pynucastro as pyna


class TestPythonNetwork2:
    @pytest.fixture(scope="class")
    def pynet(self, reaclib_library):
        mylib = reaclib_library.linking_nuclei(["he4", "c12", "o16"])
        return pyna.PythonNetwork(libraries=[mylib], inert_nuclei=["ne20"])

    def test_full_ydot_string(self, pynet):

        dyodt = \
"""dYdt[jO16] = (
   -Y[jO16]*rate_eval.O16__He4_C12
   +rho*Y[jHe4]*Y[jC12]*rate_eval.He4_C12__O16
   )

"""

        dynedt = \
"""dYdt[jNe20] = 0.0

"""

        assert pynet.full_ydot_string(pyna.Nucleus("o16")) == dyodt
        assert pynet.full_ydot_string(pyna.Nucleus("ne20")) == dynedt
