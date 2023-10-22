import filecmp
import os

import pytest

import pynucastro as pyna


@pytest.fixture(scope="package")
def reaclib_library():
    return pyna.ReacLibLibrary()


@pytest.fixture(scope="package")
def tabular_library():
    return pyna.TabularLibrary()


@pytest.fixture(scope="package")
def compare_network_files():
    def _compare_network_files(test_path, ref_path, skip_files=()):
        base_path = os.path.relpath(os.path.dirname(__file__))
        ref_path = os.path.join(base_path, ref_path)

        skip_files = set(skip_files)
        test_files = set(os.listdir(test_path)) - skip_files
        ref_files = set(os.listdir(ref_path)) - skip_files
        # files that are missing from test_path
        missing_files = ref_files - test_files
        # files that are present in test_path but not in ref_path
        extra_files = test_files - ref_files

        modified_files = set()
        # only compare files that exist in both directories, to avoid errors
        for file in test_files & ref_files:
            # note, _test is written under whatever directory pytest is run from,
            # so it is not necessarily at the same place as _amrexastro_reference
            if not filecmp.cmp(os.path.normpath(os.path.join(test_path, file)),
                               os.path.normpath(os.path.join(ref_path, file)),
                               shallow=False):
                modified_files.add(file)

        # record which files make the test fail
        if missing_files:
            print("\nmissing files:")
            for file in sorted(missing_files):
                print("  " + file)
        if extra_files:
            print("\nextra files:")
            for file in sorted(extra_files):
                print("  " + file)
        if modified_files:
            print("\nmodified files:")
            for file in sorted(modified_files):
                print("  " + file)

        assert not (
            missing_files | extra_files | modified_files
        ), "written network files don't match the stored reference"

    return _compare_network_files
