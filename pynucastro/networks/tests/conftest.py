"""A collection of fixtures for testing networks."""

import filecmp
import shutil
from pathlib import Path

import pytest

import pynucastro as pyna


@pytest.fixture(scope="package")
def reaclib_library():
    """Return the ReacLibLibrary"""
    return pyna.ReacLibLibrary()


@pytest.fixture(scope="package")
def tabular_library():
    """Return the TabularLibrary"""
    return pyna.TabularLibrary()


@pytest.fixture(scope="package")
def suzuki_library():
    """Return the Suzukibrary"""
    return pyna.SuzukiLibrary()


@pytest.fixture(scope="package")
def langanke_library():
    """Return the LangankebLibrary"""
    return pyna.LangankeLibrary()


@pytest.fixture(scope="package")
def ffn_library():
    """Return the FFNLibrary"""
    return pyna.FFNLibrary()


@pytest.fixture(scope="package")
def oda_library():
    """Return the OdaLibrary"""
    return pyna.OdaLibrary()


@pytest.fixture(scope="package")
def pruetfuller_library():
    """Return the PruetFullerLibrary"""
    return pyna.PruetFullerLibrary()


@pytest.fixture(scope="package")
def compare_network_files(request):
    """Compare files created during testing to the reference versions."""
    # this fixture returns a closure so we don't have to get the pytest config
    # in each test function and pass it through
    update_networks = request.config.getoption("--update-networks")

    def _compare_network_files(test_path, ref_path, skip_files=()):
        base_path = Path(__file__).parent.relative_to(Path.cwd())
        test_path = Path(test_path)
        ref_path = base_path/ref_path

        skip_files = set(skip_files)
        test_files = set(file.name for file in test_path.iterdir()) - skip_files
        ref_files = set(file.name for file in ref_path.iterdir()) - skip_files
        # files that are missing from test_path
        missing_files = ref_files - test_files
        # files that are present in test_path but not in ref_path
        extra_files = test_files - ref_files

        modified_files = set()
        # only compare files that exist in both directories, to avoid errors
        for file in test_files & ref_files:
            # note, _test is written under whatever directory pytest is run from,
            # so it is not necessarily at the same place as _amrexastro_reference
            if not filecmp.cmp(test_path/file, ref_path/file, shallow=False):
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
                (ref_path/file).resolve().unlink()
            # copy new and modified files to ref_path
            for file in extra_files | modified_files:
                shutil.copy(test_path/file, ref_path/file)
        else:
            assert not errors, "written network files don't match the stored reference"

    return _compare_network_files
