"""Support for a pure C++ reaction network.  These functions will
write the C++ code necessary to integrate a reaction network
comprised of the rates that are passed in.

"""


import os
import re
import shutil
import sys
from abc import ABC, abstractmethod

import numpy as np
import sympy

from pynucastro.networks.rate_collection import RateCollection
from pynucastro.networks.sympy_network_support import SympyRates


class BaseCxxNetwork(ABC, RateCollection):
    """Interpret the collection of rates and nuclei and produce the
    C++ code needed to integrate the network.

    """

    def __init__(self, *args, **kwargs):
        """Initialize the C++ network.  We take a single argument: a list
        of rate files that will make up the network

        """

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

        # a dictionary of functions to call to handle specific parts
        # of the C++ template
        self.ftags = {}
        self.ftags['<nrat_reaclib>'] = self._nrat_reaclib
        self.ftags['<nrat_tabular>'] = self._nrat_tabular
        self.ftags['<nrxn>'] = self._nrxn
        self.ftags['<rate_names>'] = self._rate_names
        self.ftags['<ebind>'] = self._ebind
        self.ftags['<compute_screening_factors>'] = self._compute_screening_factors
        self.ftags['<table_num>'] = self._table_num
        self.ftags['<declare_tables>'] = self._declare_tables
        self.ftags['<table_declare_meta>'] = self._table_declare_meta
        self.ftags['<table_init_meta>'] = self._table_init_meta
        self.ftags['<table_term_meta>'] = self._table_term_meta
        self.ftags['<compute_tabular_rates>'] = self._compute_tabular_rates
        self.ftags['<ydot>'] = self._ydot
        self.ftags['<enuc_add_energy_rate>'] = self._enuc_add_energy_rate
        self.ftags['<jacnuc>'] = self._jacnuc
        self.ftags['<initial_mass_fractions>'] = self._initial_mass_fractions
        self.ftags['<pynucastro_home>'] = self._pynucastro_home
        self.ftags['<reaclib_rate_functions>'] = self._reaclib_rate_functions
        self.ftags['<fill_reaclib_rates>'] = self._fill_reaclib_rates
        self.ftags['<approx_rate_functions>'] = self._approx_rate_functions
        self.ftags['<fill_approx_rates>'] = self._fill_approx_rates
        self.ftags['<part_fun_data>'] = self._fill_parition_function_data
        self.ftags['<part_fun_cases>'] = self._fill_parition_function_cases
        self.ftags['<spin_state_cases>'] = self._fill_spin_state_cases
        self.indent = '    '

        self.num_screen_calls = None

    @abstractmethod
    def _get_template_files(self):
        # This method should be overridden by derived classes
        # to support specific output templates.
        # This method returns a list of strings that are file paths to template files.
        return []

    def get_indent_amt(self, l, k):
        """determine the amount of spaces to indent a line"""
        rem = re.match(r'\A'+k+r'\(([0-9]*)\)\Z', l)
        return int(rem.group(1))

    def _write_network(self, odir=None):
        """
        This writes the RHS, jacobian and ancillary files for the system of ODEs that
        this network describes, using the template files.
        """
        # pylint: disable=arguments-differ

        # Prepare RHS terms
        if not self.solved_ydot:
            self.compose_ydot()
        if not self.solved_jacobian:
            self.compose_jacobian()

        # Process template files
        for tfile in self.template_files:
            tfile_basename = os.path.basename(tfile)
            outfile = tfile_basename.replace('.template', '')
            if odir is not None:
                if not os.path.isdir(odir):
                    try:
                        os.mkdir(odir)
                    except OSError:
                        sys.exit(f"unable to create directory {odir}")
                outfile = os.path.normpath(odir + "/" + outfile)

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
            tdir = os.path.dirname(tr.rfile_path)
            if tdir != os.getcwd():
                tdat_file = os.path.join(tdir, tr.table_file)
                if os.path.isfile(tdat_file):
                    shutil.copy(tdat_file, odir or os.getcwd())
                else:
                    print(f'WARNING: Table data file {tr.table_file} not found.')

    def compose_ydot(self):
        """create the expressions for dYdt for the nuclei, where Y is the
        molar fraction.


        This will take the form of a dict, where the key is a nucleus, and the
        value is a list of tuples, with the forward-reverse pairs of a rate
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
        """Create the Jacobian matrix, df/dY"""
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
        screening_map = self.get_screening_map()
        for i, scr in enumerate(screening_map):

            nuc1_info = f'{float(scr.n1.Z)}_rt, {float(scr.n1.A)}_rt'
            nuc2_info = f'{float(scr.n2.Z)}_rt, {float(scr.n2.A)}_rt'

            if not (scr.n1.dummy or scr.n2.dummy):
                # Scope the screening calculation to avoid multiple definitions of scn_fac.
                of.write(f'\n{self.indent*n_indent}' + '{')

                of.write(f'\n{self.indent*(n_indent+1)}constexpr auto scn_fac = scrn::calculate_screen_factor({nuc1_info}, {nuc2_info});\n\n')

                # Insert a static assert (which will always pass) to require the
                # compiler to evaluate the screen factor at compile time.
                of.write(f'\n{self.indent*(n_indent+1)}static_assert(scn_fac.z1 == {float(scr.n1.Z)}_rt);\n\n')

                of.write(f'\n{self.indent*(n_indent+1)}actual_screen<do_T_derivatives>(pstate, scn_fac, scor, dscor_dt);\n')

                of.write(f'{self.indent*n_indent}' + '}\n\n')

            if scr.name == "he4_he4_he4":
                # we don't need to do anything here, but we want to avoid immediately applying the screening
                pass

            elif scr.name == "he4_he4_he4_dummy":
                # make sure the previous iteration was the first part of 3-alpha
                assert screening_map[i - 1].name == "he4_he4_he4"
                # handle the second part of the screening for 3-alpha
                of.write(f'\n{self.indent*n_indent}' + '{')

                of.write(f'\n{self.indent*(n_indent+1)}constexpr auto scn_fac2 = scrn::calculate_screen_factor({nuc1_info}, {nuc2_info});\n\n')

                of.write(f'\n{self.indent*(n_indent+1)}static_assert(scn_fac2.z1 == {float(scr.n1.Z)}_rt);\n\n')

                of.write(f'\n{self.indent*(n_indent+1)}actual_screen<do_T_derivatives>(pstate, scn_fac2, scor2, dscor2_dt);\n')

                of.write(f'\n{self.indent*n_indent}' + '}\n\n')

                # there might be both the forward and reverse 3-alpha
                # if we are doing symmetric screening

                for rr in scr.rates:
                    of.write('\n')
                    of.write(f'{self.indent*n_indent}ratraw = rate_eval.screened_rates(k_{rr.fname});\n')
                    of.write(f'{self.indent*n_indent}rate_eval.screened_rates(k_{rr.fname}) *= scor * scor2;\n')
                    of.write(f'{self.indent*n_indent}if constexpr (std::is_same<T, rate_derivs_t>::value) {{\n')
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
                    of.write(f'{self.indent*n_indent}if constexpr (std::is_same<T, rate_derivs_t>::value) {{\n')
                    of.write(f'{self.indent*n_indent}    dratraw_dT = rate_eval.dscreened_rates_dT(k_{rr.fname});\n')
                    of.write(f'{self.indent*n_indent}    rate_eval.dscreened_rates_dT(k_{rr.fname}) = ratraw * dscor_dt + dratraw_dT * scor;\n')
                    of.write(f'{self.indent*n_indent}}}\n')

            of.write('\n')

        # the C++ screen.H code requires that there be at least 1 screening
        # factor because it statically allocates some arrays, so if we turned
        # off screening, just set num_screen_calls = 1 here.

        self.num_screen_calls = max(1, len(screening_map))

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

    def _table_num(self, n_indent, of):
        of.write(f'{self.indent*n_indent}const int num_tables = {len(self.tabular_rates)};\n')

    def _declare_tables(self, n_indent, of):
        for r in self.tabular_rates:
            idnt = self.indent*n_indent

            of.write(f'{idnt}extern AMREX_GPU_MANAGED table_t {r.table_index_name}_meta;\n')
            of.write(f'{idnt}extern AMREX_GPU_MANAGED Array3D<Real, 1, {r.table_temp_lines}, 1, {r.table_rhoy_lines}, 1, {r.table_num_vars}> {r.table_index_name}_data;\n')
            of.write(f'{idnt}extern AMREX_GPU_MANAGED Array1D<Real, 1, {r.table_rhoy_lines}> {r.table_index_name}_rhoy;\n')
            of.write(f'{idnt}extern AMREX_GPU_MANAGED Array1D<Real, 1, {r.table_temp_lines}> {r.table_index_name}_temp;\n')
            of.write('\n')

    def _table_declare_meta(self, n_indent, of):
        for r in self.tabular_rates:
            idnt = self.indent*n_indent

            of.write(f"{idnt}AMREX_GPU_MANAGED table_t {r.table_index_name}_meta;\n")

            of.write(f'{idnt}AMREX_GPU_MANAGED Array3D<Real, 1, {r.table_temp_lines}, 1, {r.table_rhoy_lines}, 1, {r.table_num_vars}> {r.table_index_name}_data;\n')

            of.write(f'{idnt}AMREX_GPU_MANAGED Array1D<Real, 1, {r.table_rhoy_lines}> {r.table_index_name}_rhoy;\n')
            of.write(f'{idnt}AMREX_GPU_MANAGED Array1D<Real, 1, {r.table_temp_lines}> {r.table_index_name}_temp;\n\n')

    def _table_init_meta(self, n_indent, of):
        for r in self.tabular_rates:
            idnt = self.indent*n_indent
            of.write(f'{idnt}{r.table_index_name}_meta.ntemp = {r.table_temp_lines};\n')
            of.write(f'{idnt}{r.table_index_name}_meta.nrhoy = {r.table_rhoy_lines};\n')
            of.write(f'{idnt}{r.table_index_name}_meta.nvars = {r.table_num_vars};\n')
            of.write(f'{idnt}{r.table_index_name}_meta.nheader = {r.table_header_lines};\n\n')

            of.write(f'{idnt}init_tab_info({r.table_index_name}_meta, "{r.table_file}", {r.table_index_name}_rhoy, {r.table_index_name}_temp, {r.table_index_name}_data);\n\n')

            of.write('\n')

    def _table_term_meta(self, n_indent, of):
        for r in self.tabular_rates:

            of.write('{}deallocate(num_temp_{})\n'.format(
                self.indent*n_indent, r.table_index_name))

            of.write('{}deallocate(num_rhoy_{})\n'.format(
                self.indent*n_indent, r.table_index_name))

            of.write('{}deallocate(num_vars_{})\n'.format(
                self.indent*n_indent, r.table_index_name))

            of.write('{}deallocate(rate_table_{})\n'.format(
                self.indent*n_indent, r.table_index_name))

            of.write('{}deallocate(rhoy_table_{})\n'.format(
                self.indent*n_indent, r.table_index_name))

            of.write('{}deallocate(temp_table_{})\n'.format(
                self.indent*n_indent, r.table_index_name))

            of.write('\n')

    def _compute_tabular_rates(self, n_indent, of):
        if len(self.tabular_rates) > 0:

            idnt = self.indent*n_indent

            for r in self.tabular_rates:

                of.write(f'{idnt}tabular_evaluate({r.table_index_name}_meta, {r.table_index_name}_rhoy, {r.table_index_name}_temp, {r.table_index_name}_data,\n')
                of.write(f'{idnt}                 rhoy, state.T, rate, drate_dt, edot_nu);\n')

                of.write(f'{idnt}rate_eval.screened_rates(k_{r.fname}) = rate;\n')

                of.write(f'{idnt}if constexpr (std::is_same<T, rate_derivs_t>::value) {{\n')
                of.write(f'{idnt}    rate_eval.dscreened_rates_dT(k_{r.fname}) = drate_dt;\n')
                of.write(f'{idnt}}}\n')

                of.write(f'{idnt}rate_eval.add_energy_rate(k_{r.fname}) = edot_nu;\n')
                of.write('\n')

    def _ydot(self, n_indent, of):
        # Write YDOT
        for n in self.unique_nuclei:
            if self.ydot_out_result[n] is None:
                of.write(f"{self.indent*n_indent}{self.symbol_rates.name_ydot_nuc}({n.cindex()}) = 0.0;\n\n")
                continue

            of.write(f"{self.indent*n_indent}{self.symbol_rates.name_ydot_nuc}({n.cindex()}) =\n")
            for j, pair in enumerate(self.ydot_out_result[n]):
                # pair here is the forward, reverse pair for a single rate as it affects
                # nucleus n

                if pair.count(None) == 0:
                    num = 2
                elif pair.count(None) == 1:
                    num = 1
                else:
                    raise NotImplementedError("a rate pair must contain atleast one rate")

                of.write(f"{2*self.indent*n_indent}")
                if num == 2:
                    of.write("(")

                if pair[0] is not None:
                    sol_value = self.symbol_rates.cxxify(sympy.cxxcode(pair[0], precision=15,
                                                                       standard="c++11"))

                    of.write(f"{sol_value}")

                if num == 2:
                    of.write(" + ")

                if pair[1] is not None:
                    sol_value = self.symbol_rates.cxxify(sympy.cxxcode(pair[1], precision=15,
                                                                       standard="c++11"))

                    of.write(f"{sol_value}")

                if num == 2:
                    of.write(")")

                if j == len(self.ydot_out_result[n])-1:
                    of.write(";\n\n")
                else:
                    of.write(" +\n")

    def _enuc_add_energy_rate(self, n_indent, of):
        # Add tabular per-reaction neutrino energy generation rates to the energy generation rate
        # (not thermal neutrinos)

        idnt = self.indent * n_indent

        for r in self.tabular_rates:
            if len(r.reactants) != 1:
                sys.exit('ERROR: Unknown energy rate corrections for a reaction where the number of reactants is not 1.')
            else:
                reactant = r.reactants[0]
                of.write(f'{idnt}enuc += C::Legacy::n_A * {self.symbol_rates.name_y}({reactant.cindex()}) * rate_eval.add_energy_rate(k_{r.fname});\n')

    def _jacnuc(self, n_indent, of):
        # now make the Jacobian
        n_unique_nuclei = len(self.unique_nuclei)
        for jnj, nj in enumerate(self.unique_nuclei):
            for ini, ni in enumerate(self.unique_nuclei):
                jac_idx = n_unique_nuclei*jnj + ini
                if not self.jac_null_entries[jac_idx]:
                    jvalue = self.symbol_rates.cxxify(sympy.cxxcode(self.jac_out_result[jac_idx], precision=15,
                                                                     standard="c++11"))
                    of.write(f"{self.indent*(n_indent)}scratch = {jvalue};\n")
                    of.write(f"{self.indent*n_indent}jac.set({nj.cindex()}, {ni.cindex()}, scratch);\n\n")

    def _initial_mass_fractions(self, n_indent, of):
        for i, _ in enumerate(self.unique_nuclei):
            if i == 0:
                of.write(f"{self.indent*n_indent}unit_test.X{i+1} = 1.0\n")
            else:
                of.write(f"{self.indent*n_indent}unit_test.X{i+1} = 0.0\n")

    def _pynucastro_home(self, n_indent, of):
        of.write('{}PYNUCASTRO_HOME := {}\n'.format(self.indent*n_indent,
                                                    os.path.dirname(self.pynucastro_dir)))

    def _reaclib_rate_functions(self, n_indent, of):
        assert n_indent == 0, "function definitions must be at top level"
        for r in self.reaclib_rates + self.derived_rates:
            of.write(r.function_string_cxx(dtype=self.dtype, specifiers=self.function_specifier))

    def _approx_rate_functions(self, n_indent, of):
        assert n_indent == 0, "function definitions must be at top level"
        for r in self.approx_rates:
            of.write(r.function_string_cxx(dtype=self.dtype, specifiers=self.function_specifier))

    def _fill_reaclib_rates(self, n_indent, of):
        for r in self.reaclib_rates + self.derived_rates:
            of.write(f"{self.indent*n_indent}rate_{r.fname}<do_T_derivatives>(tfactors, rate, drate_dT);\n")
            of.write(f"{self.indent*n_indent}rate_eval.screened_rates(k_{r.fname}) = rate;\n")
            of.write(f"{self.indent*n_indent}if constexpr (std::is_same<T, rate_derivs_t>::value) {{\n")
            of.write(f"{self.indent*n_indent}    rate_eval.dscreened_rates_dT(k_{r.fname}) = drate_dT;\n\n")
            of.write(f"{self.indent*n_indent}}}\n")

    def _fill_approx_rates(self, n_indent, of):
        for r in self.approx_rates:
            of.write(f"{self.indent*n_indent}rate_{r.fname}<T>(rate_eval, rate, drate_dT);\n")
            of.write(f"{self.indent*n_indent}rate_eval.screened_rates(k_{r.fname}) = rate;\n")
            of.write(f"{self.indent*n_indent}if constexpr (std::is_same<T, rate_derivs_t>::value) {{\n")
            of.write(f"{self.indent*n_indent}    rate_eval.dscreened_rates_dT(k_{r.fname}) = drate_dT;\n\n")
            of.write(f"{self.indent*n_indent}}}\n")

    def _fill_parition_function_data(self, n_indent, of):

        nuclei_pfs = self.get_nuclei_needing_partition_functions()

        decl = "MICROPHYSICS_UNUSED HIP_CONSTEXPR static AMREX_GPU_MANAGED amrex::Real"

        if nuclei_pfs:
            for n in nuclei_pfs:
                if n.partition_function:
                    npts = len(n.partition_function.temperature)
                    assert len(n.partition_function.temperature) == len(n.partition_function.partition_function)

                    # write the data out, but for readability, split it to 5 values per line

                    # number of points

                    of.write(f"{self.indent*n_indent}constexpr int {n}_npts = {npts};\n\n")

                    # first the temperature (in units of T9)

                    of.write(f"{self.indent*n_indent}{decl} {n}_temp_array[{n}_npts] = {{\n")

                    tdata = list(n.partition_function.temperature/1.0e9)
                    while tdata:
                        if len(tdata) >= 5:
                            tmp = ", ".join([str(tdata.pop(0)) for _ in range(0, 5)])
                            if tdata:
                                tmp += ","
                        else:
                            tmp = ", ".join([str(tdata.pop(0)) for _ in range(0, len(tdata))])

                        of.write(f"{2*self.indent*n_indent}{tmp}\n")
                    of.write(f"{self.indent*n_indent}}};\n\n")

                    # now the partition function itself -- this is log10

                    of.write(f"{self.indent*n_indent}// this is log10(partition function)\n\n")

                    of.write(f"{self.indent*n_indent}{decl} {n}_pf_array[{n}_npts] = {{\n")

                    tdata = list(n.partition_function.partition_function)
                    while tdata:
                        if len(tdata) >= 5:
                            tmp = ", ".join([str(np.log10(tdata.pop(0))) for _ in range(0, 5)])
                            if tdata:
                                tmp += ","
                        else:
                            tmp = ", ".join([str(np.log10(tdata.pop(0))) for _ in range(0, len(tdata))])

                        of.write(f"{2*self.indent*n_indent}{tmp}\n")
                    of.write(f"{self.indent*n_indent}}};\n\n")

                    of.write("\n")

    def _fill_parition_function_cases(self, n_indent, of):

        nuclei_pfs = self.get_nuclei_needing_partition_functions()

        if nuclei_pfs:
            for n in nuclei_pfs:
                if n.partition_function:

                    of.write(f"{self.indent*n_indent}case {n.cindex()}:\n")
                    of.write(f"{self.indent*2*n_indent}part_fun::interpolate_pf(tfactors.T9, part_fun::{n}_npts, part_fun::{n}_temp_array, part_fun::{n}_pf_array, pf, dpf_dT);\n")
                    of.write(f"{self.indent*2*n_indent}break;\n\n")

    def _fill_spin_state_cases(self, n_indent, of):

        for n in self.unique_nuclei:
            if n.spin_states is not None:
                of.write(f"{self.indent*n_indent}case {n.cindex()}:\n")
                of.write(f"{self.indent*2*n_indent}spin = {n.spin_states};\n")
                of.write(f"{self.indent*2*n_indent}break;\n\n")
        for n in self.approx_nuclei:
            if n.spin_states is not None:
                of.write(f"{self.indent*n_indent}case {n.cindex()}:\n")
                of.write(f"{self.indent*2*n_indent}spin = {n.spin_states};\n")
                of.write(f"{self.indent*2*n_indent}break;\n\n")
