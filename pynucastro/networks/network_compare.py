"""Helper functions for comparing the output of different backends
(python, C++) for the same network.  This will build and run each net
and return the ydots.
"""

import importlib
import re
import subprocess
from pathlib import Path

import numpy as np

from pynucastro.networks.amrexastro_cxx_network import AmrexAstroCxxNetwork
from pynucastro.networks.python_network import PythonNetwork
from pynucastro.networks.simple_cxx_network import SimpleCxxNetwork
from pynucastro.nucdata import Composition, Nucleus
from pynucastro.screening import chugunov_2007


class NetworkCompare:
    """A simple class to manage building a network with different
    backends to facilitate comparisons.

    Parameters
    ----------
    lib : Library
        pynucastro Library containing the rates
    use_screening : bool
        do we include screening in the comparison (Chugunov 2007)
    include_amrex : bool
        do we build and evaluate an AmrexAstroCxxNetwork?
    include_simple_cxx : bool
        do we build and evaluate a SimpleCxxNetwork?
    python_module_name : str
        the name of the module file to output the python net to
    amrex_test_path : pathlib.Path
        the Path to output the AMReX C++ source files and build in
    cxx_test_path : pathlib.Path
        the Path to output the C++ source files and build in
    """

    def __init__(self, lib, *,
                 use_screening=False,
                 include_amrex=False,
                 include_simple_cxx=True,
                 python_module_name="compare_net.py",
                 amrex_test_path=None,
                 cxx_test_path=None):

        self.lib = lib

        if use_screening:
            self.screen_func = chugunov_2007
        else:
            self.screen_func = None

        self.include_amrex = include_amrex
        self.include_simple_cxx = include_simple_cxx

        if not python_module_name.endswith(".py"):
            raise ValueError("invalid python_module_name")
        self.python_module_name = python_module_name

        if self.include_amrex:
            if amrex_test_path is None:
                amrex_test_path = Path("_test_amrex_compare/")

            self.amrex_test_path = amrex_test_path

            # make the directory where the C++ output will go
            self.amrex_test_path.mkdir(parents=True, exist_ok=True)

        if self.include_simple_cxx:
            if cxx_test_path is None:
                cxx_test_path = Path("_test_compare/")

            self.cxx_test_path = cxx_test_path

            # make the directory where the C++ output will go
            self.cxx_test_path.mkdir(parents=True, exist_ok=True)

        # we will always have the python version
        self.pynet = PythonNetwork(libraries=[lib])

        # always use a uniform compostion
        self.comp = Composition(self.pynet.unique_nuclei)
        self.comp.set_equal()

        # storage for the ydots
        self.ydots_py_inline = None
        self.ydots_py_module = None
        self.ydots_amrex = None
        self.ydots_cxx = None

        self.T_eval = None
        self.rho_eval = None

    def _run_python_inline_version(self, rho=2.e8, T=1.e9):
        """Evaluate the rates using the methods built into
        RateCollection

        """

        self.ydots_py_inline = self.pynet.evaluate_ydots(rho=rho, T=T,
                                                         composition=self.comp,
                                                         screen_func=self.screen_func)

    def _run_python_module_version(self, rho=2.e8, T=1.e9):
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
        _tmp = cn.rhs(0.0, Y, rho, T, screen_func=self.screen_func)
        self.ydots_py_module = {}
        for n, y in zip(self.pynet.unique_nuclei, _tmp):
            self.ydots_py_module[n] = y

    def _run_amrex_version(self, rho=2.e8, T=1.e9):
        """Output the AMReX C++ network code, build it, run, and
        parse the output to get the rates.

        """

        cxx_net = AmrexAstroCxxNetwork(libraries=[self.lib])
        cxx_net.write_network(odir=self.amrex_test_path, standalone_build=True)

        # build an run the simple C++ network
        if self.screen_func is not None:
            opts = "SCREEN_METHOD=chugunov2007"
        else:
            opts = "SCREEN_METHOD=null"

        # to reduce compilation time, we'll build in debug mode
        subprocess.run(f"make DEBUG=TRUE -j 4 {opts}",
                       capture_output=True,
                       shell=True, check=True,
                       cwd=self.amrex_test_path)

        cp = subprocess.run(f"./main3d.gnu.DEBUG.ex testing.density={rho} testing.temperature={T}",
                            capture_output=True,
                            shell=True, check=True, text=True,
                            cwd=self.amrex_test_path)
        stdout = cp.stdout

        # the stdout includes lines of the form:
        #    Ydot(X) = ...
        # for each nucleus X.  This regex will capture
        # the nucleus and the ydot value for each of these
        ydot_re = re.compile(r"(Ydot)\((\w*)\)(\s+)(=)(\s+)([\d\-e\+.]*)",
                             re.IGNORECASE | re.DOTALL)

        self.ydots_amrex = {}
        for line in stdout.split("\n"):
            if match := ydot_re.search(line.strip()):
                self.ydots_amrex[Nucleus(match.group(2))] = float(match.group(6))

    def _run_simple_cxx_version(self, rho=2.e8, T=1.e9):
        """Output the simple C++ network code, build it, run, and
        parse the output to get the rates.

        """

        cxx_net = SimpleCxxNetwork(libraries=[self.lib])
        cxx_net.write_network(odir=self.cxx_test_path)

        # build an run the simple C++ network
        if self.screen_func is not None:
            opts = "USE_SCREENING=TRUE"
        else:
            opts = "USE_SCREENING=FALSE"

        subprocess.run(f"make -j 4 {opts}",
                       capture_output=False,
                       shell=True, check=True,
                       cwd=self.cxx_test_path)

        cp = subprocess.run(f"./main {rho} {T}",
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
                self.ydots_cxx[Nucleus(match.group(2))] = float(match.group(6))

    def evaluate(self, rho=2.e8, T=1.e9):
        """Evaluate the ydots from all the backends we are
        considering

        Parameters
        ----------
        rho : float
            density to evaluate rates at (CGS)
        T : float
            temperature to evaluate rates at (K)

        """

        self._run_python_inline_version(rho=rho, T=T)
        self._run_python_module_version(rho=rho, T=T)

        if self.include_amrex:
            self._run_amrex_version(rho=rho, T=T)

        if self.include_simple_cxx:
            self._run_simple_cxx_version(rho=rho, T=T)

        self.T_eval = T
        self.rho_eval = rho

    def print_summary(self):
        """Print a summary of the dY/dt comparison and errors for each
        network type run.

        """

        # we need to have run previously
        if self.ydots_py_inline is None:
            raise ValueError("no ydots stored.  evaluate() must be run first")

        data_headers = {"py (inline)": self.ydots_py_inline,
                        "py (module)": self.ydots_py_module}

        if self.ydots_amrex:
            data_headers["AMReX C++"] = self.ydots_amrex

        if self.ydots_cxx:
            data_headers["simple C++"] = self.ydots_cxx

        header = f" {'nuc':5} "
        for key in data_headers:
            if key == "py (inline)":
                header += f"| {key:13} "
            else:
                header += f"| {key:13} {'error':11} "

        print(header)
        print("-" * len(header))

        for nuc in self.ydots_py_inline:
            line = f" {nuc!s:5} "
            for key, ydots in data_headers.items():
                val = ydots[nuc]
                ref = self.ydots_py_inline[nuc]
                if key == "py (inline)":
                    line += f"| {val:13.6g} "
                else:
                    err = abs((val - ref) / ref)
                    line += f"| {val:13.6g} {err:11.5g} "
            print(line)
