import shutil
from pathlib import Path

import numpy as np
import pytest
from numpy.testing import assert_allclose

from pynucastro import networks


class TestFullPythonNetwork:
    @pytest.fixture(scope="class")
    def fn(self, reaclib_library, suzuki_library):
        rate_names = ["c12(c12,a)ne20",
                      "c12(c12,n)mg23",
                      "c12(c12,p)na23",
                      "c12(a,g)o16",
                      "n(,)p",
                      "he4(aa,g)c12"]
        rates = reaclib_library.get_rate_by_name(rate_names)

        tabular_rate_names = ["na23(,)ne23",
                              "ne23(,)na23"]
        tabular_rates = suzuki_library.get_rate_by_name(tabular_rate_names)

        fn = networks.PythonNetwork(rates=rates+tabular_rates)
        return fn

    def test_write_network(self, fn, compare_network_files):
        """test the write_network function"""
        test_path = Path("_test_python/")
        # subdirectory of pynucastro/networks/tests/
        reference_path = Path("_python_reference/")
        # files that will be ignored if present in the generated directory
        skip_files = []

        test_file = "network.py"

        # remove any previously generated files
        shutil.rmtree(test_path, ignore_errors=True)
        test_path.mkdir(parents=True, exist_ok=True)
        fn.write_network(outfile=test_path/test_file)
        compare_network_files(test_path, reference_path, skip_files)

        # clean up generated files if the test passed
        shutil.rmtree(test_path)

    def test_ydots(self, fn):

        test_path = Path("_test_python/")
        test_file = "network.py"

        test_path.mkdir(parents=True, exist_ok=True)
        fn.write_network(outfile=test_path/test_file)

        import _test_python.network as net

        X = np.zeros(net.nnuc)
        X[:] = 1.0 / net.nnuc
        Y = X * net.nnuc

        rho = 1.e8
        T = 1.e9

        ydot = net.rhs(0.0, Y, rho, T)

        ydot_benchmark = np.array([-1.129734e-03,  1.979093e-03, -1.702699e+06,  5.667057e+05,
                                   6.454310e+02,  1.044892e-03, -1.1838042e-2,  1.268268e-02,
                                   4.712856e-06])

        assert_allclose(ydot, ydot_benchmark, rtol=1.e-6)
