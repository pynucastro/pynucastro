import filecmp
import os
import shutil
import warnings

import pytest

import pynucastro as pyna


@pytest.fixture(scope="package")
def reaclib_library():
    return pyna.ReacLibLibrary()


@pytest.fixture(scope="package")
def tabular_library():
    return pyna.TabularLibrary()


@pytest.fixture(scope="package")
def compare_network_files(request):
    # this fixture returns a closure so we don't have to get the pytest config
    # in each test function and pass it through
    update_networks = request.config.getoption("--update-networks")

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

        errors = missing_files | extra_files | modified_files
        if update_networks and errors:
            # remove files that are no longer present in test_path
            for file in missing_files:
                os.unlink(os.path.normpath(os.path.join(ref_path, file)))
            # copy new and modified files to ref_path
            for file in extra_files | modified_files:
                shutil.copy(os.path.normpath(os.path.join(test_path, file)),
                            os.path.normpath(os.path.join(ref_path, file)))
            # raise a warning to let the user know which tests were updated
            warnings.warn(UserWarning(f"updated reference files in {ref_path}"))
        else:
            assert not errors, "written network files don't match the stored reference"

    return _compare_network_files
