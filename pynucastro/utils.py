"""Various utilities for working with pynucastro."""


import json
import subprocess
from importlib.metadata import Distribution
from urllib.parse import urlparse

from ._version import version


def pynucastro_version():
    """Get as descriptive of a version as possible.  If this was
    installed via pip directly, then it will be the version assigned
    to that release.  However, if we are installed in editable mode,
    then the version from __version__ might not be updated, since that
    is created dynamically at the pip invocation.  In that case, we
    will use git to find the version.

    Returns
    -------
    str

    """

    direct_url = Distribution.from_name("pynucastro").read_text("direct_url.json")
    if direct_url is None:
        # in pytest, things are run in a new environment without a
        # distribution
        return version

    is_editable = json.loads(direct_url).get("dir_info", {}).get("editable", False)

    if is_editable:
        # we need to use git to find the full version
        git_dir = json.loads(direct_url).get("url", None)
        if git_dir:
            git_dir = urlparse(git_dir).path

            try:
                sp = subprocess.run("git describe", capture_output=True,
                                    shell=True, check=True, text=True,
                                    cwd=git_dir)
                git_version = sp.stdout.strip()
            except subprocess.CalledProcessError:
                git_version = version

            return git_version

    return version
