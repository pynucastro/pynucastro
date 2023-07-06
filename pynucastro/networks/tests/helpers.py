import filecmp
import os


def compare_network_files(test_path, ref_path, skip_files=()):
    # get the list of files from the generated directory, so any new files
    # will need to be added to the reference directory or explicitly ignored
    files = os.listdir(test_path)

    errors = []
    for test_file in files:
        if test_file in skip_files:
            continue
        # note, _test is written under whatever directory pytest is run from,
        # so it is not necessarily at the same place as _amrexastro_reference
        if not filecmp.cmp(os.path.normpath(f"{test_path}/{test_file}"),
                           os.path.normpath(f"{ref_path}/{test_file}"),
                           shallow=False):
            errors.append(test_file)

    return errors
