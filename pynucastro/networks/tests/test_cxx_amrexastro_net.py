# unit tests for rates
import io
import shutil

import pytest

from pynucastro import networks


class TestAmrexAstroCxxNetwork:
    # pylint: disable=protected-access
    @pytest.fixture(scope="class")
    def fn(self, reaclib_library, suzuki_library):
        rate_names = ["c12(c12,a)ne20",
                      "c12(c12,n)mg23",
                      "c12(c12,p)na23",
                      "c12(a,g)o16",
                      "n(,)p"]
        rates = reaclib_library.get_rate_by_name(rate_names)

        tabular_rate_names = ["na23(,)ne23",
                              "ne23(,)na23"]
        tabular_rates = suzuki_library.get_rate_by_name(tabular_rate_names)

        fn = networks.AmrexAstroCxxNetwork(rates=rates+tabular_rates)
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

        answer = ('    k_C12_C12_to_He4_Ne20 = 1,\n' +
                  '    k_C12_C12_to_n_Mg23 = 2,\n' +
                  '    k_C12_C12_to_p_Na23 = 3,\n' +
                  '    k_He4_C12_to_O16 = 4,\n' +
                  '    k_n_to_p_weak_wc12 = 5,\n' +
                  '    k_Na23_to_Ne23 = 6,\n' +
                  '    k_Ne23_to_Na23 = 7,\n' +
                  '    NumRates = k_Ne23_to_Na23\n')
        assert self.cromulent_ftag(fn._nrxn, answer, n_indent=1)

    def test_ebind(self, fn):
        """ test the _ebind function """

        answer = ('        if constexpr (spec == N) {\n' +
                  '            return 0.0_rt;\n' +
                  '        }\n' +
                  '        else if constexpr (spec == H1) {\n' +
                  '            return 0.0_rt;\n' +
                  '        }\n' +
                  '        else if constexpr (spec == He4) {\n' +
                  '            return 28.29566_rt;\n' +
                  '        }\n' +
                  '        else if constexpr (spec == C12) {\n' +
                  '            return 92.16172800000001_rt;\n' +
                  '        }\n' +
                  '        else if constexpr (spec == O16) {\n' +
                  '            return 127.619296_rt;\n' +
                  '        }\n' +
                  '        else if constexpr (spec == Ne20) {\n' +
                  '            return 160.6448_rt;\n' +
                  '        }\n' +
                  '        else if constexpr (spec == Ne23) {\n' +
                  '            return 182.970888_rt;\n' +
                  '        }\n' +
                  '        else if constexpr (spec == Na23) {\n' +
                  '            return 186.56433900000002_rt;\n' +
                  '        }\n' +
                  '        else if constexpr (spec == Mg23) {\n' +
                  '            return 181.725645_rt;\n' +
                  '        }\n')

        assert self.cromulent_ftag(fn._ebind, answer, n_indent=2)

    def test_write_network(self, fn, compare_network_files):
        """ test the write_network function"""
        test_path = "_test_cxx/"
        # subdirectory of pynucastro/networks/tests/
        reference_path = "_amrexastro_cxx_reference/"
        # files that will be ignored if present in the generated directory
        skip_files = []

        # remove any previously generated files
        shutil.rmtree(test_path, ignore_errors=True)
        fn.write_network(odir=test_path)
        compare_network_files(test_path, reference_path, skip_files)

        # clean up generated files if the test passed
        shutil.rmtree(test_path)
