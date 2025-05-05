# this creates the same network as SimpleCxxNetwork and PythonNetwork
# and have them both compute dY/dt and compares to make sure that they
# agree.  Note: screening is not considered.

import os
import re
import shlex
import subprocess
from pathlib import Path

import numpy as np
import pytest
from pytest import approx

import pynucastro as pyna


def run(string):
    prog = shlex.split(string)
    p0 = subprocess.Popen(prog, stdin=None, stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE)
    stdout0, stderr0 = p0.communicate()
    rc = p0.returncode
    p0.stdout.close()
    p0.stderr.close()
    return rc, stdout0.decode("utf-8"), stderr0.decode("utf-8")


class TestNetworkCompare:

    @pytest.fixture(scope="class")
    def lib(self, reaclib_library):
        nuc = ["p", "he4", "c12", "o16", "ne20", "na23", "mg24"]
        lib = reaclib_library.linking_nuclei(nuc)
        return lib

    def test_compare(self, lib):

        test_path = Path("_test_compare/")

        test_path.mkdir(parents=True, exist_ok=True)
        os.chdir(test_path)

        # C++

        cxx_net = pyna.SimpleCxxNetwork(libraries=[lib])
        cxx_net.write_network()

        rc, _, _ = run("make DISABLE_SCREENING=TRUE")
        assert rc == 0

        _, stdout, _ = run("./main")

        ydot_re = re.compile(r"(Ydot)\((\w*)\)(\s+)(=)(\s+)([\d\-e\+.]*)",
                             re.IGNORECASE | re.DOTALL)

        ydots_cxx = {}
        for line in stdout.split("\n"):
            if match := ydot_re.search(line.strip()):
                ydots_cxx[match.group(2)] = float(match.group(6))

        os.chdir("../")

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
