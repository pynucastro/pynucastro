import filecmp
import os


def compare_network_files(test_path, ref_path, skip_files=()):
    base_path = os.path.relpath(os.path.dirname(__file__))
    ref_path = os.path.join(base_path, ref_path)

    skip_files = set(skip_files)
    test_files = set(os.listdir(test_path)) - skip_files
    ref_files = set(os.listdir(ref_path)) - skip_files
    assert test_files == ref_files, "missing/extra files"

    errors = []
    for file in sorted(test_files):
        # note, _test is written under whatever directory pytest is run from,
        # so it is not necessarily at the same place as _amrexastro_reference
        if not filecmp.cmp(os.path.normpath(os.path.join(test_path, file)),
                           os.path.normpath(os.path.join(ref_path, file)),
                           shallow=False):
            errors.append(file)

    assert not errors, f"files don't match: {' '.join(errors)}"
