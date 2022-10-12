# unit tests for rates
import filecmp
import io
import os

import pytest

import pynucastro as pyna


class TestAmrexAstroCxxNetwork:
    # pylint: disable=protected-access
    @pytest.fixture(scope="class")
    def fn(self, reaclib_library):

        mylib = reaclib_library.linking_nuclei(["mg24", "al27", "si28", "p31", "s32", "he4", "p"])
        net = pyna.AmrexAstroCxxNetwork(libraries=[mylib])
        net.make_ap_pg_approx()
        net.remove_nuclei(["al27", "p31"])
        fn = net
        return fn

    def cromulent_ftag(self, ftag, answer, n_indent=1):
        """ check to see if function ftag returns answer """

        output = io.StringIO()
        ftag(n_indent, output)
        result = output.getvalue() == answer
        output.close()
        return result

    def test_write_network(self, fn):
        """ test the write_network function"""
        test_path = "_test_cxx_approx/"
        reference_path = "_amrexastro_cxx_approx_reference/"
        base_path = os.path.relpath(os.path.dirname(__file__))

        fn.write_network(odir=test_path)

        files = ["actual_network_data.cpp",
                 "actual_network.H",
                 "actual_rhs.H",
                 "inputs.burn_cell.VODE",
                 "Make.package",
                 "NETWORK_PROPERTIES",
                 "_parameters",
                 "pynucastro.net",
                 "reaclib_rates.H",
                 "table_rates_data.cpp",
                 "table_rates.H"]

        errors = []
        for test_file in files:
            # note, _test is written under whatever directory pytest is run from,
            # so it is not necessarily at the same place as _amrexastro_reference
            if not filecmp.cmp(os.path.normpath(f"{test_path}/{test_file}"),
                               os.path.normpath(f"{base_path}/{reference_path}/{test_file}"),
                               shallow=False):
                errors.append(test_file)

        assert not errors, f"files don't match: {' '.join(errors)}"
