# unit tests for rates
import shutil

import pytest

import pynucastro as pyna
from pynucastro.rates.alternate_rates import IliadisO16pgF17


class TestAmrexAstroCxxNetwork:
    @pytest.fixture(scope="class")
    def fn(self):

        iliadis = IliadisO16pgF17()
        iliadis_derived = pyna.DerivedRate(iliadis, use_pf=True)
        fn = pyna.AmrexAstroCxxNetwork(rates=[iliadis, iliadis_derived])
        return fn

    def test_write_network(self, fn, compare_network_files):
        """ test the write_network function"""
        test_path = "_test_cxx_temptab/"
        # subdirectory of pynucastro/networks/tests/
        reference_path = "_amrexastro_cxx_temptab_reference/"
        # files that will be ignored if present in the generated directory
        skip_files = ["pynucastro-info.txt"]

        # remove any previously generated files
        shutil.rmtree(test_path, ignore_errors=True)
        fn.write_network(odir=test_path)
        compare_network_files(test_path, reference_path, skip_files)

        # clean up generated files if the test passed
        shutil.rmtree(test_path)
