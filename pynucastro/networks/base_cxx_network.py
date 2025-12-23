"""Support for a pure C++ reaction network.  These functions will
write the C++ code necessary to integrate a reaction network
comprised of the rates that are passed in.

"""


import itertools
import re
import shutil
import sys
import warnings
from abc import ABC, abstractmethod
from pathlib import Path

import numpy as np
import sympy

from pynucastro.constants import constants
from pynucastro.networks.rate_collection import RateCollection
from pynucastro.networks.sympy_network_support import SympyRates
from pynucastro.screening import get_screening_map
from pynucastro.utils import pynucastro_version


def _rate_dtype(nrxn):
    """Given the number of reactions (nrxn), return the smallest C++
    unsigned integer type that can hold them

    """

    dtype = "std::uint32_t"
    # give 1 extra padding in case we use a final value in an enum
    if nrxn < 255:
        dtype = "std::uint8_t"
    elif nrxn < 65535:
        dtype = "std::uint16_t"
    return dtype


def _signed_rate_dtype(nrxn):
    """Given the number of reactions (nrxn), return the smallest C++
    signed integer type that can hold them

    """

    dtype = "int"
    if nrxn < 127:
        dtype = "std::int8_t"
    elif nrxn < 32767:
        dtype = "short"
    return dtype


class BaseCxxNetwork(ABC, RateCollection):
    """Base class for a C++ network.  This takes the same arguments as
    :py:class:`RateCollection
    <pynucastro.networks.rate_collection.RateCollection>` and
    interprets the collection of rates and nuclei to produce the C++
    code needed to integrate the network.

    """

    def __init__(self, *args, **kwargs):

        super().__init__(*args, **kwargs)

        # Get the template files for writing this network code
        self.template_files = self._get_template_files()

        self.symbol_rates = SympyRates()

        self.ydot_out_result = None
        self.solved_ydot = False
        self.jac_out_result = None
        self.jac_null_entries = None
        self.solved_jacobian = False

        self.function_specifier = "inline"
        self.dtype = "double"
        self.array_namespace = ""

        # a dictionary of functions to call to handle specific parts
        # of the C++ template
        self.ftags = {}
        self.ftags['<nrat_reaclib>'] = self._nrat_reaclib
        self.ftags['<nrat_tabular>'] = self._nrat_tabular
        self.ftags['<nrxn>'] = self._nrxn
        self.ftags['<nrxn_enum_type>'] = self._nrxn_enum_type
        self.ftags['<rate_names>'] = self._rate_names
        self.ftags['<ebind>'] = self._ebind
        self.ftags['<mion>'] = self._mion
        self.ftags['<compute_screening_factors>'] = self._compute_screening_factors
        self.ftags['<table_num>'] = self._table_num
        self.ftags['<declare_tables>'] = self._declare_tables
        self.ftags['<table_declare_meta>'] = self._table_declare_meta
        self.ftags['<table_init_meta>'] = self._table_init_meta
        self.ftags['<compute_tabular_rates>'] = self._compute_tabular_rates
        self.ftags['<ydot>'] = self._ydot
        self.ftags['<ydot_weak>'] = self._ydot_weak
        self.ftags['<jacnuc>'] = self._jacnuc
        self.ftags['<reaclib_rate_functions>'] = self._reaclib_rate_functions
        self.ftags['<rate_struct>'] = self._rate_struct
        self.ftags['<fill_reaclib_rates>'] = self._fill_reaclib_rates
        self.ftags['<derived_rate_functions>'] = self._derived_rate_functions
        self.ftags['<fill_derived_rates>'] = self._fill_derived_rates
        self.ftags['<approx_rate_functions>'] = self._approx_rate_functions
        self.ftags['<fill_approx_rates>'] = self._fill_approx_rates
        self.ftags['<part_fun_data>'] = self._fill_partition_function_data
        self.ftags['<part_fun_declare>'] = self._fill_partition_function_declare
        self.ftags['<part_fun_cases>'] = self._fill_partition_function_cases
        self.ftags['<declare_pf_cache_temp_index>'] = self._declare_pf_cache_temp_index
        self.ftags['<spin_state_cases>'] = self._fill_spin_state_cases
        self.ftags['<pynucastro_version>'] = self._fill_pynucastro_version
        self.indent = '    '

    @abstractmethod
    def _get_template_files(self):
        # This method should be overridden by derived classes
        # to support specific output templates.
        # This method returns a list of strings that are file paths to template files.
        return []

    def get_indent_amt(self, l, k):
        """Determine the amount of spaces to indent a line.  This
        looks for a string of the form "<key>(#)", where # is the a
        number that is the amount of indent.

        Parameters
        ----------
        l : str
            a line from a template file to check for an indent
        k : str
            a keyword of the form "<key>" that appears in the line

        Returns
        -------
        int

        """

        rem = re.match(r'\A'+k+r'\(([0-9]*)\)\Z', l)
        return int(rem.group(1))

    def _write_network(self, odir=None):
        """Output the RHS, jacobian and ancillary files for the system
        of ODEs that this network describes, using the template files.

        """
        # pylint: disable=arguments-differ

        # Prepare RHS terms
        if not self.solved_ydot:
            self.compose_ydot()
        if not self.solved_jacobian:
            self.compose_jacobian()

        # Process template files
        for tfile in self.template_files:
            outfile = tfile.name.replace('.template', '')
            if odir is not None:
                odir = Path(odir)
                if not odir.is_dir():
                    try:
                        odir.mkdir()
                    except OSError:
                        sys.exit(f"unable to create directory {odir}")
                outfile = odir/outfile

            with open(tfile) as ifile, open(outfile, "w") as of:
                for l in ifile:
                    ls = l.strip()
                    foundkey = False
                    for k, func in self.ftags.items():
                        if k in ls:
                            foundkey = True
                            n_indent = self.get_indent_amt(ls, k)
                            func(n_indent, of)
                    if not foundkey:
                        of.write(l)

        # Copy any tables in the network to the current directory
        # if the table file cannot be found, print a warning and continue.
        for tr in self.tabular_rates:
            tdir = tr.rfile_path.resolve().parent
            if tdir != Path.cwd():
                tdat_file = Path(tdir, tr.table_file)
                if tdat_file.is_file():
                    shutil.copy(tdat_file, odir or Path.cwd())
                else:
                    warnings.warn(UserWarning(f'Table data file {tr.table_file} not found.'))

    def compose_ydot(self):
        """Create the expressions for dY/dt for each nucleus, where Y
        is the molar fraction.

        This stores the result in a dict where the key is a
        :py:class:`Nucleus <pynucastro.nucdata.nucleus.Nucleus>`, and
        the value is a list of tuples, with the forward-reverse pairs
        of a rate

        """

        ydot = {}
        for n in self.unique_nuclei:
            if not self.nuclei_rate_pairs[n]:
                ydot[n] = None
            else:
                ydot_sym_terms = []
                for rp in self.nuclei_rate_pairs[n]:
                    if rp.forward is not None:
                        fwd = self.symbol_rates.ydot_term_symbol(rp.forward, n)
                    else:
                        fwd = None

                    if rp.reverse is not None:
                        rvs = self.symbol_rates.ydot_term_symbol(rp.reverse, n)
                    else:
                        rvs = None

                    ydot_sym_terms.append((fwd, rvs))
                ydot[n] = ydot_sym_terms

        self.ydot_out_result = ydot
        self.solved_ydot = True

    def compose_jacobian(self):
        """Create the Jacobian matrix, df/dY, where f is a dY/dt and Y
        is a molar fraction

        The Jacobian is stored as a list with each entry representing
        a Jacobian element.  We also store whether the entry is null.

        """
        jac_null = []
        jac_sym = []
        for nj in self.unique_nuclei:
            for ni in self.unique_nuclei:
                rsym_is_null = True
                rsym = float(sympy.sympify(0.0))
                for r in self.nuclei_consumed[nj]:
                    rsym_add, rsym_add_null = self.symbol_rates.jacobian_term_symbol(r, nj, ni)
                    rsym = rsym + rsym_add
                    rsym_is_null = rsym_is_null and rsym_add_null
                for r in self.nuclei_produced[nj]:
                    rsym_add, rsym_add_null = self.symbol_rates.jacobian_term_symbol(r, nj, ni)
                    rsym = rsym + rsym_add
                    rsym_is_null = rsym_is_null and rsym_add_null
                jac_sym.append(rsym)
                jac_null.append(rsym_is_null)

        self.jac_out_result = jac_sym
        self.jac_null_entries = jac_null
        self.solved_jacobian = True

    def _compute_screening_factors(self, n_indent, of):
        if not self.do_screening:
            screening_map = []
        else:
            screening_map = get_screening_map(self.get_rates())
        for i, scr in enumerate(screening_map):

            nuc1_info = f'{float(scr.n1.Z)}_rt, {float(scr.n1.A)}_rt'
            nuc2_info = f'{float(scr.n2.Z)}_rt, {float(scr.n2.A)}_rt'

            if not (scr.n1.dummy or scr.n2.dummy):
                # Scope the screening calculation to avoid multiple definitions of scn_fac.
                of.write(f'\n{self.indent*n_indent}' + '{\n')

                of.write(f'{self.indent*(n_indent+1)}constexpr auto scn_fac = scrn::calculate_screen_factor({nuc1_info}, {nuc2_info});\n')

                # Insert a static assert (which will always pass) to require the
                # compiler to evaluate the screen factor at compile time.
                of.write(f'{self.indent*(n_indent+1)}static_assert(scn_fac.z1 == {float(scr.n1.Z)}_rt);\n')

                of.write(f'{self.indent*(n_indent+1)}actual_screen(pstate, scn_fac, scor, dscor_dt);\n')

                of.write(f'{self.indent*n_indent}' + '}\n')

            if scr.name == "He4_He4_He4":
                # we don't need to do anything here, but we want to avoid immediately applying the screening
                pass

            elif scr.name == "He4_He4_He4_dummy":
                # make sure the previous iteration was the first part of 3-alpha
                assert screening_map[i - 1].name == "He4_He4_He4"
                # handle the second part of the screening for 3-alpha
                of.write(f'\n{self.indent*n_indent}' + '{\n')

                of.write(f'{self.indent*(n_indent+1)}constexpr auto scn_fac2 = scrn::calculate_screen_factor({nuc1_info}, {nuc2_info});\n')

                of.write(f'{self.indent*(n_indent+1)}static_assert(scn_fac2.z1 == {float(scr.n1.Z)}_rt);\n')

                of.write(f'{self.indent*(n_indent+1)}actual_screen(pstate, scn_fac2, scor2, dscor2_dt);\n')

                of.write(f'{self.indent*n_indent}' + '}\n')

                # we can have both a(aa,g)c12 and a(aa,p)b11
                for rr in scr.rates:
                    of.write('\n')
                    of.write(f'{self.indent*n_indent}ratraw = rate_eval.screened_rates(k_{rr.fname});\n')
                    of.write(f'{self.indent*n_indent}rate_eval.screened_rates(k_{rr.fname}) *= scor * scor2;\n')
                    of.write(f'{self.indent*n_indent}if constexpr (std::is_same_v<T, rate_derivs_t>) {{\n')
                    of.write(f'{self.indent*n_indent}    dratraw_dT = rate_eval.dscreened_rates_dT(k_{rr.fname});\n')
                    of.write(f'{self.indent*n_indent}    rate_eval.dscreened_rates_dT(k_{rr.fname}) = ratraw * (scor * dscor2_dt + dscor_dt * scor2) + dratraw_dT * scor * scor2;\n')
                    of.write(f'{self.indent*n_indent}}}\n')
            else:
                # there might be several rates that have the same
                # reactants and therefore the same screening applies
                # -- handle them all now

                for rr in scr.rates:
                    of.write('\n')
                    of.write(f'{self.indent*n_indent}ratraw = rate_eval.screened_rates(k_{rr.fname});\n')
                    of.write(f'{self.indent*n_indent}rate_eval.screened_rates(k_{rr.fname}) *= scor;\n')
                    of.write(f'{self.indent*n_indent}if constexpr (std::is_same_v<T, rate_derivs_t>) {{\n')
                    of.write(f'{self.indent*n_indent}    dratraw_dT = rate_eval.dscreened_rates_dT(k_{rr.fname});\n')
                    of.write(f'{self.indent*n_indent}    rate_eval.dscreened_rates_dT(k_{rr.fname}) = ratraw * dscor_dt + dratraw_dT * scor;\n')
                    of.write(f'{self.indent*n_indent}}}\n')

            of.write('\n')

    def _nrat_reaclib(self, n_indent, of):
        # Writes the number of Reaclib rates
        of.write(f'{self.indent*n_indent}const int NrateReaclib = {len(self.reaclib_rates + self.derived_rates)};\n')

    def _nrat_tabular(self, n_indent, of):
        # Writes the number of tabular rates
        of.write(f'{self.indent*n_indent}const int NrateTabular = {len(self.tabular_rates)};\n')

    def _nrxn(self, n_indent, of):
        for i, r in enumerate(self.all_rates):
            of.write(f'{self.indent*n_indent}k_{r.fname} = {i+1},\n')
        of.write(f'{self.indent*n_indent}NumRates = k_{self.all_rates[-1].fname}\n')

    def _nrxn_enum_type(self, n_indent, of):
        nrxn = len(self.all_rates)
        dtype = _rate_dtype(nrxn)
        of.write(f'{self.indent*n_indent}{dtype}\n')

    def _rate_names(self, n_indent, of):
        for i, r in enumerate(self.all_rates):
            if i < len(self.all_rates)-1:
                cont = ","
            else:
                cont = ""
            of.write(f'{self.indent*n_indent}"{r.fname}"{cont}  // {i+1},\n')

    def _ebind(self, n_indent, of):
        for nuc in self.unique_nuclei:
            of.write(f'{self.indent*n_indent}ebind_per_nucleon({nuc.cindex()}) = {nuc.nucbind}_rt;\n')

    def _mion(self, n_indent, of):
        for nuc in self.unique_nuclei:
            of.write(f'{self.indent*n_indent}mion({nuc.cindex()}) = {nuc.A_nuc * constants.m_u_C18}_rt;\n')

    def _table_num(self, n_indent, of):
        of.write(f'{self.indent*n_indent}const int num_tables = {len(self.tabular_rates)};\n')

    def _declare_tables(self, n_indent, of):
        for r in self.tabular_rates:
            idnt = self.indent*n_indent

            of.write(f'{idnt}extern AMREX_GPU_MANAGED table_t {r.table_index_name}_meta;\n')
            of.write(f'{idnt}extern AMREX_GPU_MANAGED {self.array_namespace}Array3D<{self.dtype}, 1, {r.table_temp_lines}, 1, {r.table_rhoy_lines}, 1, {r.table_num_vars}> {r.table_index_name}_data;\n')
            of.write(f'{idnt}extern AMREX_GPU_MANAGED {self.array_namespace}Array1D<{self.dtype}, 1, {r.table_rhoy_lines}> {r.table_index_name}_rhoy;\n')
            of.write(f'{idnt}extern AMREX_GPU_MANAGED {self.array_namespace}Array1D<{self.dtype}, 1, {r.table_temp_lines}> {r.table_index_name}_temp;\n')
            of.write('\n')

    def _table_declare_meta(self, n_indent, of):
        for r in self.tabular_rates:
            idnt = self.indent*n_indent

            of.write(f"{idnt}AMREX_GPU_MANAGED table_t {r.table_index_name}_meta;\n")

            of.write(f'{idnt}AMREX_GPU_MANAGED {self.array_namespace}Array3D<{self.dtype}, 1, {r.table_temp_lines}, 1, {r.table_rhoy_lines}, 1, {r.table_num_vars}> {r.table_index_name}_data;\n')

            of.write(f'{idnt}AMREX_GPU_MANAGED {self.array_namespace}Array1D<{self.dtype}, 1, {r.table_rhoy_lines}> {r.table_index_name}_rhoy;\n')
            of.write(f'{idnt}AMREX_GPU_MANAGED {self.array_namespace}Array1D<{self.dtype}, 1, {r.table_temp_lines}> {r.table_index_name}_temp;\n\n')

    def _table_init_meta(self, n_indent, of):
        for r in self.tabular_rates:
            idnt = self.indent*n_indent
            of.write(f'{idnt}{r.table_index_name}_meta.ntemp = {r.table_temp_lines};\n')
            of.write(f'{idnt}{r.table_index_name}_meta.nrhoy = {r.table_rhoy_lines};\n')
            of.write(f'{idnt}{r.table_index_name}_meta.nvars = {r.table_num_vars};\n')
            of.write(f'{idnt}{r.table_index_name}_meta.nheader = {r.table_header_lines};\n\n')

            of.write(f'{idnt}init_tab_info({r.table_index_name}_meta, "{r.rfile}", {r.table_index_name}_rhoy, {r.table_index_name}_temp, {r.table_index_name}_data);\n\n')

            of.write('\n')

    def _compute_tabular_rates(self, n_indent, of):
        if len(self.tabular_rates) > 0:

            idnt = self.indent*n_indent

            of.write(f'{idnt}amrex::Real log_temp = std::log10(state.T);\n')
            of.write(f'{idnt}amrex::Real log_rhoy = std::log10(rhoy);\n\n')

            for r in self.tabular_rates:

                of.write(f'{idnt}tabular_evaluate({r.table_index_name}_meta, {r.table_index_name}_rhoy, {r.table_index_name}_temp, {r.table_index_name}_data,\n')
                of.write(f'{idnt}                 log_rhoy, log_temp, state.T, rate, drate_dt, edot_nu, edot_gamma);\n')

                of.write(f'{idnt}rate_eval.screened_rates(k_{r.fname}) = rate;\n')

                of.write(f'{idnt}if constexpr (std::is_same_v<T, rate_derivs_t>) {{\n')
                of.write(f'{idnt}    rate_eval.dscreened_rates_dT(k_{r.fname}) = drate_dt;\n')
                of.write(f'{idnt}}}\n')

                of.write(f'{idnt}rate_eval.enuc_weak += C::n_A * {self.symbol_rates.name_y}({r.reactants[0].cindex()}) * (edot_nu + edot_gamma);\n')

                of.write('\n')

    def _cxxify(self, s):
        # This is a helper function that converts sympy cxxcode to the actual c++ code we use.
        return self.symbol_rates.cxxify(s)

    def _write_ydot_nuc(self, n_indent, of, ydot_nuc):
        # Helper function to write out ydot of a specific nuclei

        for j, pair in enumerate(ydot_nuc):
            # pair here is the forward, reverse pair for a single rate as it affects
            # nucleus n

            if pair.count(None) == 0:
                num = 2
            elif pair.count(None) == 1:
                num = 1
            else:
                raise NotImplementedError("a rate pair must contain at least one rate")

            of.write(f"{2*self.indent*n_indent}")
            if num == 2:
                of.write("(")

            if pair[0] is not None:
                sol_value = self._cxxify(sympy.cxxcode(pair[0], precision=15,
                                                       standard="c++11"))

                of.write(f"{sol_value}")

            if num == 2:
                of.write(" + ")

            if pair[1] is not None:
                sol_value = self._cxxify(sympy.cxxcode(pair[1], precision=15,
                                                       standard="c++11"))

                of.write(f"{sol_value}")

            if num == 2:
                of.write(")")

            if j == len(ydot_nuc)-1:
                of.write(";\n\n")
            else:
                of.write(" +\n")

    def _ydot(self, n_indent, of):
        # Write YDOT
        for n in self.unique_nuclei:
            if self.ydot_out_result[n] is None:
                of.write(f"{self.indent*n_indent}{self.symbol_rates.name_ydot_nuc}({n.cindex()}) = 0.0_rt;\n\n")
                continue

            of.write(f"{self.indent*n_indent}{self.symbol_rates.name_ydot_nuc}({n.cindex()}) =\n")

            self._write_ydot_nuc(n_indent, of, self.ydot_out_result[n])

    def _ydot_weak(self, n_indent, of):
        # Writes ydot for tabular weak reactions only

        # Get the tabular weak rates first.
        idnt = self.indent*n_indent

        if len(self.tabular_rates) > 0:

            of.write(f'{idnt}amrex::Real log_temp = std::log10(state.T);\n')
            of.write(f'{idnt}amrex::Real log_rhoy = std::log10(rhoy);\n\n')

            for r in self.tabular_rates:

                of.write(f'{idnt}tabular_evaluate({r.table_index_name}_meta, {r.table_index_name}_rhoy, {r.table_index_name}_temp, {r.table_index_name}_data,\n')
                of.write(f'{idnt}                 log_rhoy, log_temp, state.T, rate, drate_dt, edot_nu, edot_gamma);\n')

                of.write(f'{idnt}rate_eval.screened_rates(k_{r.fname}) = rate;\n')

                of.write(f'{idnt}rate_eval.enuc_weak += C::n_A * {self.symbol_rates.name_y}({r.reactants[0].cindex()}) * (edot_nu + edot_gamma);\n')

                of.write('\n')
            of.write(f'{idnt}auto screened_rates = rate_eval.screened_rates;\n')
        of.write('\n')

        # Compose and write ydot weak

        for n in self.unique_nuclei:

            has_weak_rates = any(
                (rp.forward is not None and rp.forward.weak) or
                (rp.reverse is not None and rp.reverse.weak)
                for rp in self.nuclei_rate_pairs[n]
            )

            if not self.nuclei_rate_pairs[n] or not has_weak_rates:
                of.write(f"{self.indent*n_indent}{self.symbol_rates.name_ydot_nuc}({n.cindex()}) = 0.0_rt;\n\n")
                continue

            ydot_sym_terms = []
            for rp in self.nuclei_rate_pairs[n]:
                fwd = None
                if rp.forward is not None and rp.forward.weak:
                    fwd = self.symbol_rates.ydot_term_symbol(rp.forward, n)

                rvs = None
                if rp.reverse is not None and rp.reverse.weak:
                    rvs = self.symbol_rates.ydot_term_symbol(rp.reverse, n)

                if (fwd, rvs).count(None) < 2:
                    ydot_sym_terms.append((fwd, rvs))

            of.write(f"{self.indent*n_indent}{self.symbol_rates.name_ydot_nuc}({n.cindex()}) =\n")

            self._write_ydot_nuc(n_indent, of, ydot_sym_terms)

    def _jacnuc(self, n_indent, of):
        # now make the Jacobian
        n_unique_nuclei = len(self.unique_nuclei)
        for jnj, nj in enumerate(self.unique_nuclei):
            for ini, ni in enumerate(self.unique_nuclei):
                jac_idx = n_unique_nuclei*jnj + ini
                if not self.jac_null_entries[jac_idx]:
                    jvalue = self._cxxify(sympy.cxxcode(self.jac_out_result[jac_idx], precision=15,
                                                                     standard="c++11"))
                    of.write(f"{self.indent*(n_indent)}scratch = {jvalue};\n")
                    of.write(f"{self.indent*n_indent}jac.set({nj.cindex()}, {ni.cindex()}, scratch);\n\n")
                else:
                    of.write(f"{self.indent*n_indent}jac.set({nj.cindex()}, {ni.cindex()}, 0.0);\n\n")

    def _reaclib_rate_functions(self, n_indent, of):
        assert n_indent == 0, "function definitions must be at top level"
        for r in self.reaclib_rates + self.modified_rates:
            of.write(r.function_string_cxx(dtype=self.dtype, specifiers=self.function_specifier))

    def _derived_rate_functions(self, n_indent, of):
        assert n_indent == 0, "function definitions must be at top level"
        for r in self.derived_rates:
            of.write(r.function_string_cxx(dtype=self.dtype, specifiers=self.function_specifier))

    def _rate_struct(self, n_indent, of):
        assert n_indent == 0, "function definitions must be at top level"

        of.write("struct rate_t {\n")
        of.write(f"    {self.array_namespace}Array1D<{self.dtype}, 1, NumRates>  screened_rates;\n")
        of.write(f"    {self.dtype} enuc_weak;\n")
        of.write("};\n\n")
        of.write("struct rate_derivs_t {\n")
        of.write(f"    {self.array_namespace}Array1D<{self.dtype}, 1, NumRates>  screened_rates;\n")
        of.write(f"    {self.array_namespace}Array1D<{self.dtype}, 1, NumRates>  dscreened_rates_dT;\n")
        of.write(f"    {self.dtype} enuc_weak;\n")
        of.write("};\n\n")

    def _approx_rate_functions(self, n_indent, of):
        assert n_indent == 0, "function definitions must be at top level"
        for r in self.approx_rates:
            of.write(r.function_string_cxx(dtype=self.dtype, specifiers=self.function_specifier))

    def _fill_reaclib_rates(self, n_indent, of):
        # note: modified_rates needs to be on the end here, since they
        # likely will call the underlying reaclib rate for the actual
        # rate evaluation
        for r in self.reaclib_rates + self.modified_rates:
            of.write(f"{self.indent*n_indent}rate_{r.fname}<do_T_derivatives>(tfactors, rate, drate_dT);\n")
            of.write(f"{self.indent*n_indent}rate_eval.screened_rates(k_{r.fname}) = rate;\n")
            of.write(f"{self.indent*n_indent}if constexpr (std::is_same_v<T, rate_derivs_t>) {{\n")
            of.write(f"{self.indent*n_indent}    rate_eval.dscreened_rates_dT(k_{r.fname}) = drate_dT;\n\n")
            of.write(f"{self.indent*n_indent}}}\n")

    def _fill_derived_rates(self, n_indent, of):
        if self.derived_rates:
            of.write(f"{self.indent*n_indent}part_fun::pf_cache_t pf_cache{{}};\n\n")
            temp_arrays, _ = self.dedupe_partition_function_temperatures()
            for i in range(len(temp_arrays)):
                of.write(f"{self.indent*n_indent}pf_cache.index_temp_array_{i+1} = interp_net::find_index(tfactors.T9, part_fun::temp_array_{i+1});\n")
                of.write("\n")

        for r in self.derived_rates:
            of.write(f"{self.indent*n_indent}rate_{r.fname}<T>(rate_eval, tfactors, rate, drate_dT, pf_cache);\n")
            of.write(f"{self.indent*n_indent}rate_eval.screened_rates(k_{r.fname}) = rate;\n")
            of.write(f"{self.indent*n_indent}if constexpr (std::is_same_v<T, rate_derivs_t>) {{\n")
            of.write(f"{self.indent*n_indent}    rate_eval.dscreened_rates_dT(k_{r.fname}) = drate_dT;\n\n")
            of.write(f"{self.indent*n_indent}}}\n")

    def _fill_approx_rates(self, n_indent, of):
        for r in self.approx_rates:
            args = ["rate_eval"]
            if r.rate_eval_needs_rho:
                args.append("rho")
            if r.rate_eval_needs_comp:
                args.append("Y")
            args += ["rate", "drate_dT"]

            of.write(f"{self.indent*n_indent}rate_{r.fname}<T>({', '.join(args)});\n")
            of.write(f"{self.indent*n_indent}rate_eval.screened_rates(k_{r.fname}) = rate;\n")
            of.write(f"{self.indent*n_indent}if constexpr (std::is_same_v<T, rate_derivs_t>) {{\n")
            of.write(f"{self.indent*n_indent}    rate_eval.dscreened_rates_dT(k_{r.fname}) = drate_dT;\n\n")
            of.write(f"{self.indent*n_indent}}}\n")

    def _fill_partition_function_declare(self, n_indent, of):

        temp_arrays, temp_indices = self.dedupe_partition_function_temperatures()

        for i, temp in enumerate(temp_arrays):

            decl = f"extern AMREX_GPU_MANAGED amrex::Array1D<{self.dtype}, 0, npts_{i+1}-1>"

            # number of points
            of.write(f"{self.indent*n_indent}constexpr int npts_{i+1} = {len(temp)};\n\n")

            # write the temperature array sizes out

            of.write(f"{self.indent*n_indent}// this is T9\n\n")

            of.write(f"{self.indent*n_indent}{decl} temp_array_{i+1};\n\n")

        for n, i in temp_indices.items():
            # declare the partition function data

            of.write(f"{self.indent*n_indent}// this is log10(partition function)\n\n")

            decl = f"extern AMREX_GPU_MANAGED amrex::Array1D<{self.dtype}, 0, npts_{i+1}-1>"
            of.write(f"{self.indent*n_indent}{decl} {n}_pf_array;\n")
            thresh_temp = n.get_part_func_threshold_temp()
            # convert to T9 if it is physical
            if thresh_temp > 0:
                thresh_temp /= 1.e9
            of.write(f"{self.indent*n_indent}constexpr {self.dtype} {n}_pf_threshold_T9 = {thresh_temp};\n\n")

    def _fill_partition_function_data(self, n_indent, of):
        # itertools recipe
        def batched(iterable, n):
            """Batch data into tuples of length n. The last batch may
            be shorter.

            """
            # batched('ABCDEFG', 3) --> ABC DEF G
            if n < 1:
                raise ValueError('n must be at least one')
            it = iter(iterable)
            while batch := tuple(itertools.islice(it, n)):
                yield batch

        temp_arrays, temp_indices = self.dedupe_partition_function_temperatures()

        for i, temp in enumerate(temp_arrays):
            # number of points

            decl = f"AMREX_GPU_MANAGED amrex::Array1D<{self.dtype}, 0, npts_{i+1}-1>"

            # write the temperature out, but for readability, split it to 5 values per line

            of.write(f"{self.indent*n_indent}// this is T9\n\n")

            of.write(f"{self.indent*n_indent}{decl} temp_array_{i+1}= {{\n")

            for data in batched(temp / 1.0e9, 5):
                tmp = " ".join([f"{t}," for t in data])
                of.write(f"{self.indent*(n_indent+1)}{tmp}\n")
            of.write(f"{self.indent*n_indent}}};\n\n")

            if i == len(temp_arrays) - 1:
                of.write("\n")

        for n, i in temp_indices.items():
            # write the partition function data out, but for readability, split
            # it to 5 values per line

            of.write(f"{self.indent*n_indent}// this is log10(partition function)\n\n")

            decl = f"AMREX_GPU_MANAGED amrex::Array1D<{self.dtype}, 0, npts_{i+1}-1>"
            of.write(f"{self.indent*n_indent}{decl} {n}_pf_array = {{\n")

            for data in batched(np.log10(n.partition_function.partition_function), 5):
                tmp = " ".join([f"{x}," for x in data])
                of.write(f"{self.indent*(n_indent+1)}{tmp}\n")
            of.write(f"{self.indent*n_indent}}};\n\n")

    def _fill_partition_function_cases(self, n_indent, of):

        _, temp_indices = self.dedupe_partition_function_temperatures()

        for n, i in temp_indices.items():
            of.write(f"{self.indent*n_indent}case {n.cindex()}:\n")
            of.write(f"{self.indent*(n_indent+1)}if (tfactors.T9 > part_fun::{n}_pf_threshold_T9) {{\n")
            of.write(f"{self.indent*(n_indent+2)}part_fun::interpolate_pf(tfactors.T9, pf_cache.index_temp_array_{i+1}, part_fun::temp_array_{i+1}, part_fun::{n}_pf_array, pf, dpf_dT);\n")
            of.write(f"{self.indent*(n_indent+1)}}}\n")
            of.write(f"{self.indent*(n_indent+1)}break;\n\n")

    def _declare_pf_cache_temp_index(self, n_indent, of):

        temp_arrays, _ = self.dedupe_partition_function_temperatures()

        for i, _ in enumerate(temp_arrays):

            # number of points
            of.write(f"{self.indent*n_indent}int index_temp_array_{i+1}{{-1}};\n\n")

    def _fill_spin_state_cases(self, n_indent, of):

        def key_func(nuc):
            if nuc.spin_states is None:
                return -1
            return nuc.spin_states

        # group identical cases together to satisfy clang-tidy
        nuclei = sorted(self.unique_nuclei + self.approx_nuclei, key=key_func)
        for spin_state, group in itertools.groupby(nuclei, key=key_func):
            if spin_state == -1:
                continue
            for n in group:
                of.write(f"{self.indent*n_indent}case {n.cindex()}:\n")
            of.write(f"{self.indent*(n_indent+1)}spin = {spin_state};\n")
            of.write(f"{self.indent*(n_indent+1)}break;\n\n")

    def _fill_pynucastro_version(self, n_indent, of):
        of.write(f"{self.indent*n_indent}pynucastro version: {pynucastro_version()}\n")
