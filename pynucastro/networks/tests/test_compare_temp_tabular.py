# this creates a network with a TemperatureTabularRate and compares to
# ensure that both the RateCollection version and the PythonNetwork
# written to a module give the same ydots

import numpy as np
import pytest
from pytest import approx

import pynucastro as pyna
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

        # python
        pynet = pyna.PythonNetwork(libraries=[lib])
        pynet.write_network("temptab_net.py")

        import temptab_net as tt  # pylint: disable=import-outside-toplevel,import-error

        rho = 1.e5
        T = 2.e8
        comp = pyna.Composition(pynet.unique_nuclei)
        comp.set_equal()

        # compare to the RateCollection version

        ydots_py = pynet.evaluate_ydots(rho=rho, T=T, composition=comp)

        Y = np.asarray(list(comp.get_molar().values()))
        module_ydots = tt.rhs(0.0, Y, rho, T)
        for n, k in enumerate(ydots_py):
            assert ydots_py[k] == approx(module_ydots[n], rel=1.e-11, abs=1.e-14)
