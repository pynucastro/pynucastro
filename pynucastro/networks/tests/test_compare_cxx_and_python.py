# this creates the same network as SimpleCxxNetwork and PythonNetwork.
# we then compare the ydots from the C++ net, the python network
# written to a module, and the python network evaluating as a
# RateCollection.  Note: screening is not considered.

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

        # thermodynamic conditions
        # we set the composition to be uniform for all tests
        rho = 2.e8
        T = 1.e9

        # METHOD 1:
        # python / RateCollection version -- computed via
        # pynucastro eval methods.  This serves as the baseline

        pynet = pyna.PythonNetwork(libraries=[lib])

        comp = pyna.Composition(pynet.unique_nuclei)
        comp.set_equal()

        ydots_py = pynet.evaluate_ydots(rho=rho, T=T, composition=comp)

        # METHOD 2:
        # simple-C++ network

        cxx_net = pyna.SimpleCxxNetwork(libraries=[lib])
        cxx_net.write_network(odir=test_path)

        # build an run the simple C++ network
        subprocess.run("make USE_SCREENING=FALSE",
                       capture_output=False,
                       shell=True, check=True, cwd=test_path)

        cp = subprocess.run(f"./main {rho} {T}",
                            capture_output=True,
                            shell=True, check=True, text=True, cwd=test_path)
        stdout = cp.stdout

        # the stdout includes lines of the form:
        #    Ydot(X) = ...
        # for each nucleus X.  This regex will capture
        # the nucleus and the ydot value for each of these
        ydot_re = re.compile(r"(Ydot)\((\w*)\)(\s+)(=)(\s+)([\d\-e\+.]*)",
                             re.IGNORECASE | re.DOTALL)

        ydots_cxx = {}
        for line in stdout.split("\n"):
            if match := ydot_re.search(line.strip()):
                ydots_cxx[match.group(2)] = float(match.group(6))

        # do the comparison

        for k, v in ydots_cxx.items():
            nuc = pyna.Nucleus(k)
            assert v == approx(ydots_py[nuc], rel=1.e-11, abs=1.e-14)

        # METHOD 3:
        # python module written as a file, with a function written
        # for computing each rate

        pynet.write_network("compare_net.py")

        import compare_net as cn  # pylint: disable=import-outside-toplevel,import-error

        # we can now compute the ydots via cn.rhs()
        # do the comparison

        Y = np.asarray(list(comp.get_molar().values()))
        module_ydots = cn.rhs(0.0, Y, rho, T)
        for n, k in enumerate(ydots_cxx):
            assert ydots_cxx[k] == approx(module_ydots[n], rel=1.e-11, abs=1.e-14)
