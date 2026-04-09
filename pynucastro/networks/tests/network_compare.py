"""Helper functions for comparing the output of different backends
(python, C++) for the same network.  This will build and run each net
and return the ydots.

Note: this disables screening.

"""

import importlib
import re
import subprocess
from pathlib import Path

import numpy as np

import pynucastro as pyna


class NetworkCompare:
    """A simple class to manage building a network with different
    backends to facilitate comparisons.

    Parameters
    ----------
    lib : Library
        pynucastro Library containing the rates
    rho : float
        density to evaluate rates at (CGS)
    T : float
        temperature to evaluate rates at (K)
    include_simple_cxx : bool
        do we build and evaluate a SimpleCxxNetwork?
    python_module_name : str
        the name of the module file to output the python net to
    cxx_test_path : Path
        the Path to output the C++ source files and build in
    """

    def __init__(self, lib, *,
                 rho=2.e8, T=1.e9,
                 include_simple_cxx=True,
                 python_module_name="compare_net.py",
                 cxx_test_path=None):

        self.lib = lib
        self.rho = rho
        self.T = T

        self.include_simple_cxx = include_simple_cxx

        if not python_module_name.endswith(".py"):
            raise ValueError("invalid python_module_name")
        self.python_module_name = python_module_name

        if self.include_simple_cxx:
            if cxx_test_path is None:
                cxx_test_path = Path("_test_compare/")

            self.cxx_test_path = cxx_test_path

            # make the directory where the C++ output will go
            self.cxx_test_path.mkdir(parents=True, exist_ok=True)

        # we will always have the python version
        self.pynet = pyna.PythonNetwork(libraries=[lib])

        # always use a uniform compostion
        self.comp = pyna.Composition(self.pynet.unique_nuclei)
        self.comp.set_equal()

        # storage for the ydots
        self.ydots_py_inline = None
        self.ydots_py_module = None
        self.ydots_cxx = None

    def _run_python_inline_version(self):
        """Evaluate the rates using the methods built into
        RateCollection

        """

        self.ydots_py_inline = self.pynet.evaluate_ydots(rho=self.rho, T=self.T,
                                                         composition=self.comp)

    def _run_python_module_version(self):
        """Write the python network to a module and import it, and
        then evaluate the rates from the rate functions in the module.

        """

        self.pynet.write_network(self.python_module_name)

        spec = importlib.util.spec_from_file_location("cn",
                                                      self.python_module_name)
        cn = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(cn)

        # we can now compute the ydots via cn.rhs()
        Y = np.asarray(list(self.comp.get_molar().values()))
        _tmp = cn.rhs(0.0, Y, self.rho, self.T)
        self.ydots_py_module = {}
        for n, y in zip(self.pynet.unique_nuclei, _tmp):
            self.ydots_py_module[n] = y

    def _run_simple_cxx_version(self):
        """Output the simple C++ network code, build it, run, and
        parse the output to get the rates.

        """

        cxx_net = pyna.SimpleCxxNetwork(libraries=[self.lib])
        cxx_net.write_network(odir=self.cxx_test_path)

        # build an run the simple C++ network
        subprocess.run("make USE_SCREENING=FALSE",
                       capture_output=False,
                       shell=True, check=True,
                       cwd=self.cxx_test_path)

        cp = subprocess.run(f"./main {self.rho} {self.T}",
                            capture_output=True,
                            shell=True, check=True, text=True,
                            cwd=self.cxx_test_path)
        stdout = cp.stdout

        # the stdout includes lines of the form:
        #    Ydot(X) = ...
        # for each nucleus X.  This regex will capture
        # the nucleus and the ydot value for each of these
        ydot_re = re.compile(r"(Ydot)\((\w*)\)(\s+)(=)(\s+)([\d\-e\+.]*)",
                             re.IGNORECASE | re.DOTALL)

        self.ydots_cxx = {}
        for line in stdout.split("\n"):
            if match := ydot_re.search(line.strip()):
                self.ydots_cxx[pyna.Nucleus(match.group(2))] = float(match.group(6))

    def evaluate(self):
        """Evaluate the ydots from all the backends we are
        considering

        """

        self._run_python_inline_version()
        self._run_python_module_version()

        if self.include_simple_cxx:
            self._run_simple_cxx_version()
