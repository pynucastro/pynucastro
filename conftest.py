"""
Recent versions of numpy use AVX-512 routines from SVML to speed up some math
functions on processors that support it. These implementations allow a maximum
error of 4 ULP, which is looser than in previous versions. Our tests for
write_network() (and probably some others) expect exactly reproducible floating
point values, which may differ between the AVX-512 and fallback
implementations. Unfortunately, not all of Github's runners support AVX-512, so
we get unpredictable CI failures depending on what hardware the tests run on.

For now, we're disabling these features when running the test suite with the
NPY_DISABLE_CPU_FEATURES environment variable. This needs to be set before
numpy is imported for the first time, and will emit a RuntimeWarning if the
disabled features are not supported by the machine.

This may be fixed by https://github.com/numpy/numpy/pull/24006, which uses the
high accuracy (within 1 ULP) SVML routines for float64 values.

Related numpy PRs and issues:
  https://github.com/numpy/numpy/pull/19478
  https://github.com/numpy/numpy/pull/20991
  https://github.com/numpy/numpy/pull/22240
  https://github.com/numpy/numpy/issues/23523
"""
import os
import sys
import warnings

if sys.platform == "linux" or sys.platform == "linux2":
    os.environ["NPY_DISABLE_CPU_FEATURES"] = "AVX512F AVX512CD AVX512_SKX"
elif sys.platform == "darwin":
    os.environ["NPY_DISABLE_CPU_FEATURES"] = "ASIMDHP"

# ignore all NPY_DISABLE_CPU_FEATURES warnings in any subprocesses
# need this for nbval as it expects stderr to be empty
os.environ["PYTHONWARNINGS"] = "ignore:During parsing environment variable 'NPY_DISABLE_CPU_FEATURES':RuntimeWarning::"

with warnings.catch_warnings():
    # ignore just the "not supported by your machine" warning for the standard pytest tests
    warnings.filterwarnings(
        action="ignore",
        message=r"(?s)During parsing environment variable 'NPY_DISABLE_CPU_FEATURES'.*not supported by your machine",
        category=RuntimeWarning
    )
    import numpy  # noqa[F401]  # pylint: disable=unused-import


def pytest_addoption(parser):
    parser.addoption(
        "--update-networks",
        action="store_true",
        help="Update the reference outputs for all failing write_network tests",
    )
