# this creates a network with a TemperatureTabularRate and compares to
# ensure that both the RateCollection version and the PythonNetwork
# written to a module give the same ydots

import pytest
from pytest import approx

from pynucastro.networks.network_compare import NetworkCompare
from pynucastro.rates.alternate_rates import IliadisO16pgF17


class TestNetworkCompare:

    @pytest.fixture(scope="class")
    def lib(self, reaclib_library):
        nuc = ["p", "n15", "o16"]
        lib = reaclib_library.linking_nuclei(nuc, with_reverse=False)
        r = IliadisO16pgF17()
        lib.add_rate(r)
        return lib

    def test_compare(self, lib):

        # thermodynamic conditions
        # we set the composition to be uniform for all tests
        rho = 1.e5
        T = 2.e8

        nc = NetworkCompare(lib, rho=rho, T=T,
                            include_simple_cxx=False,
                            python_module_name="net_tt.py")
        nc.evaluate()

        # compare the inline and module python versions

        for nuc in nc.ydots_py_inline:
            assert nc.ydots_py_inline[nuc] == approx(nc.ydots_py_module[nuc],
                                                     rel=1.e-11, abs=1.e-14)
