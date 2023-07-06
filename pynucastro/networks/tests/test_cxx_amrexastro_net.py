# unit tests for rates
import io
import os.path

import pytest

from pynucastro import networks
from pynucastro.networks.tests.helpers import compare_network_files


class TestAmrexAstroCxxNetwork:
    # pylint: disable=protected-access
    @pytest.fixture(scope="class")
    def fn(self):
        files = ["c12-c12a-ne20-cf88",
                 "c12-c12n-mg23-cf88",
                 "c12-c12p-na23-cf88",
                 "c12-ag-o16-nac2",
                 "na23--ne23-toki",
                 "ne23--na23-toki",
                 "n--p-wc12"]

        fn = networks.AmrexAstroCxxNetwork(files)
        return fn

    def cromulent_ftag(self, ftag, answer, n_indent=1):
        """ check to see if function ftag returns answer """

        output = io.StringIO()
        ftag(n_indent, output)
        result = output.getvalue() == answer
        output.close()
        return result

    def test_nrat_reaclib(self, fn):
        """ test the _nrat_reaclib function """

        answer = '    const int NrateReaclib = 5;\n'

        assert self.cromulent_ftag(fn._nrat_reaclib, answer, n_indent=1)

    def test_nrat_tabular(self, fn):
        """ test the _nrat_tabular function """

        answer = '    const int NrateTabular = 2;\n'

        assert self.cromulent_ftag(fn._nrat_tabular, answer, n_indent=1)

    def test_nrxn(self, fn):
        """ test the _nrxn function """

        answer = ('    k_c12_c12_to_he4_ne20 = 1,\n' +
                  '    k_c12_c12_to_n_mg23 = 2,\n' +
                  '    k_c12_c12_to_p_na23 = 3,\n' +
                  '    k_he4_c12_to_o16 = 4,\n' +
                  '    k_n_to_p_weak_wc12 = 5,\n' +
                  '    k_na23_to_ne23 = 6,\n' +
                  '    k_ne23_to_na23 = 7,\n' +
                  '    NumRates = k_ne23_to_na23\n')
        assert self.cromulent_ftag(fn._nrxn, answer, n_indent=1)

    def test_ebind(self, fn):
        """ test the _ebind function """

        answer = ('        ebind_per_nucleon(N) = 0.0_rt;\n' +
                  '        ebind_per_nucleon(H1) = 0.0_rt;\n' +
                  '        ebind_per_nucleon(He4) = 7.073915_rt;\n' +
                  '        ebind_per_nucleon(C12) = 7.680144_rt;\n' +
                  '        ebind_per_nucleon(O16) = 7.976206_rt;\n' +
                  '        ebind_per_nucleon(Ne20) = 8.03224_rt;\n' +
                  '        ebind_per_nucleon(Ne23) = 7.955256_rt;\n' +
                  '        ebind_per_nucleon(Na23) = 8.111493000000001_rt;\n' +
                  '        ebind_per_nucleon(Mg23) = 7.901115_rt;\n')
        assert self.cromulent_ftag(fn._ebind, answer, n_indent=2)

    def test_write_network(self, fn):
        """ test the write_network function"""
        test_path = "_test_cxx/"
        reference_path = "_amrexastro_cxx_reference/"
        # files that will be ignored if present in the generated directory
        skip_files = [
            "23Na-23Ne_electroncapture.dat",
            "23Ne-23Na_betadecay.dat",
        ]

        base_path = os.path.relpath(os.path.dirname(__file__))
        reference_path = os.path.join(base_path, reference_path)

        fn.write_network(odir=test_path)
        compare_network_files(test_path, reference_path, skip_files)
