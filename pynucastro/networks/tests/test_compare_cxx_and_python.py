# this creates the same network as SimpleCxxNetwork and PythonNetwork
# and have them both compute dY/dt and compares to make sure that they
# agree.  Note: screening is not considered.

import re
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest
from pytest import approx

import pynucastro as pyna


class TestNetworkCompare:

    @pytest.fixture(scope="class")
    def lib(self, reaclib_library):
        nuc = ["p", "he4", "c12", "o16", "ne20", "na23", "mg24"]
        lib = reaclib_library.linking_nuclei(nuc)
        return lib

    @pytest.mark.skipif(sys.platform.startswith("win"), reason="Does not run on Windows")
    def test_compare(self, lib):

        test_path = Path("_test_compare/")

        test_path.mkdir(parents=True, exist_ok=True)

        # C++

        cxx_net = pyna.SimpleCxxNetwork(libraries=[lib])
        cxx_net.write_network(odir=test_path)

        subprocess.run("make DISABLE_SCREENING=TRUE", capture_output=False,
                       shell=True, check=True, cwd=test_path)

        cp = subprocess.run("./main", capture_output=True,
                            shell=True, check=True, text=True, cwd=test_path)
        stdout = cp.stdout

        ydot_re = re.compile(r"(Ydot)\((\w*)\)(\s+)(=)(\s+)([\d\-e\+.]*)",
                             re.IGNORECASE | re.DOTALL)

        ydots_cxx = {}
        for line in stdout.split("\n"):
            if match := ydot_re.search(line.strip()):
                ydots_cxx[match.group(2)] = float(match.group(6))

        # python
        pynet = pyna.PythonNetwork(libraries=[lib])
        pynet.write_network("compare_net.py")

        import compare_net as cn  # pylint: disable=import-outside-toplevel,import-error

        rho = 2.e8
        T = 1.e9
        comp = pyna.Composition(pynet.unique_nuclei)
        comp.set_equal()

        # compare to the RateCollection version

        ydots_py = pynet.evaluate_ydots(rho=rho, T=T, composition=comp)

        for k, v in ydots_cxx.items():
            nuc = pyna.Nucleus(k)
            assert v == approx(ydots_py[nuc], rel=1.e-11, abs=1.e-14)

        # compare to the module version

        Y = np.asarray(list(comp.get_molar().values()))
        module_ydots = cn.rhs(0.0, Y, rho, T)
        for n, k in enumerate(ydots_cxx):
            assert ydots_cxx[k] == approx(module_ydots[n], rel=1.e-11, abs=1.e-14)
