import os

import pytest

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
