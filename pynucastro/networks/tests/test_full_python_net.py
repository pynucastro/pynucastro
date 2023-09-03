import os
import shutil

import numpy as np
import pytest
from numpy.testing import assert_allclose

from pynucastro import networks


class TestFullPythonNetwork:
    @pytest.fixture(scope="class")
    def fn(self):
        files = ["c12-c12a-ne20-cf88",
                 "c12-c12n-mg23-cf88",
                 "c12-c12p-na23-cf88",
                 "c12-ag-o16-nac2",
                 "na23--ne23-toki",
                 "ne23--na23-toki",
                 "n--p-wc12",
                 "he4-aag-c12-fy05"]

        return networks.PythonNetwork(files)

    def test_write_network(self, fn):
        """test the write_network function"""
        test_path = "_test_python/"
        reference_path = "_python_reference/"
        base_path = os.path.relpath(os.path.dirname(__file__))

        test_file = "network.py"

        os.makedirs(test_path, exist_ok=True)
        fn.write_network(outfile=os.path.join(test_path, test_file))

        # compare contents of files
        with open(os.path.join(test_path, test_file), "r") as generated, \
             open(os.path.join(base_path, reference_path, test_file), "r") as reference:
            assert generated.readlines() == reference.readlines()

        # clean up generated files if the test passed
        shutil.rmtree(test_path)

    def test_ydots(self, fn):

        test_path = "_test_python/"
        test_file = "network.py"

        os.makedirs(test_path, exist_ok=True)
        fn.write_network(outfile=os.path.join(test_path, test_file))

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
