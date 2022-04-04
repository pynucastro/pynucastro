"""Support for a pure C++ reaction network.  These functions will
write the C++ code necessary to integrate a reaction network
comprised of the rates that are passed in.

"""


import os
import shutil
import sys
import re
from collections import OrderedDict
from abc import ABC, abstractmethod
import random
import string

import sympy
from pynucastro.networks import RateCollection
from pynucastro.networks import SympyRates

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

        self.symbol_rates = SympyRates(ctype="C++")

        self.ydot_out_result  = None
        self.solved_ydot      = False
        self.jac_out_result   = None
        self.jac_null_entries = None
        self.solved_jacobian  = False

        self.secret_code = ''.join(random.choices(string.ascii_uppercase + string.digits, k=32))

        # a dictionary of functions to call to handle specific parts
        # of the C++ template
        self.ftags = OrderedDict()
        self.ftags['<nrat_reaclib>'] = self._nrat_reaclib
        self.ftags['<nrat_tabular>'] = self._nrat_tabular
        self.ftags['<nrxn>'] = self._nrxn
        self.ftags['<ebind>'] = self._ebind
        self.ftags['<screen_add>'] = self._screen_add
        self.ftags['<compute_screening_factors>'] = self._compute_screening_factors
        self.ftags['<write_reaclib_metadata>'] = self._write_reaclib_metadata
        self.ftags['<table_num>'] = self._table_num
        self.ftags['<declare_tables>'] = self._declare_tables
        self.ftags['<table_declare_meta>'] = self._table_declare_meta
        self.ftags['<table_init_meta>'] = self._table_init_meta
        self.ftags['<table_term_meta>'] = self._table_term_meta
        self.ftags['<table_rates_indices>'] = self._table_rates_indices
        self.ftags['<compute_tabular_rates>'] = self._compute_tabular_rates
        self.ftags['<ydot>'] = self._ydot
        self.ftags['<enuc_add_energy_rate>'] = self._enuc_add_energy_rate
        self.ftags['<jacnuc>'] = self._jacnuc
        self.ftags['<initial_mass_fractions>'] = self._initial_mass_fractions
        self.ftags['<pynucastro_home>'] = self._pynucastro_home
        self.ftags['<secret_code>'] = self._secret_code_write
        self.ftags['<secret_code_set>'] = self._secret_code_write_reference
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
        rem = re.match(r'\A'+k+r'\(([0-9]*)\)\Z',l)
        return int(rem.group(1))

    def _write_network(self, odir=None):
        """
        This writes the RHS, jacobian and ancillary files for the system of ODEs that
        this network describes, using the template files.
        """

        # Prepare RHS terms
        if not self.solved_ydot:
            self.compose_ydot()
        if not self.solved_jacobian:
            self.compose_jacobian()

        # Process template files
        for tfile in self.template_files:
            tfile_basename = os.path.basename(tfile)
            outfile    = tfile_basename.replace('.template', '')
            if odir is not None:
                if not os.path.isdir(odir):
                    try:
                        os.mkdir(odir)
                    except:
                        sys.exit(f"unable to create directory {odir}")
                outfile = os.path.normpath(odir + "/" + outfile)

            with open(tfile) as ifile, open(outfile, "w") as of:
                for l in ifile:
                    ls = l.strip()
                    foundkey = False
                    for k in self.ftags:
                        if k in ls:
                            foundkey = True
                            n_indent = self.get_indent_amt(ls, k)
                            self.ftags[k](n_indent, of)
                    if not foundkey:
                        of.write(l)

        # Copy any tables in the network to the current directory
        # if the table file cannot be found, print a warning and continue.
        for i_tab in self.tabular_rates:
            tr = self.rates[i_tab]
            tdir = os.path.dirname(tr.rfile_path)
            if tdir != os.getcwd():
                tdat_file = os.path.join(tdir, tr.table_file)
                if os.path.isfile(tdat_file):
                    shutil.copy(tdat_file, os.getcwd())
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

        self.ydot_out_result  = ydot
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

        self.jac_out_result  = jac_sym
        self.jac_null_entries = jac_null
        self.solved_jacobian = True

    def _compute_screening_factors(self, n_indent, of):
        screening_map = self.get_screening_map()
        for i, scr in enumerate(screening_map):

            nuc1_info = f'{float(scr.n1.Z)}_rt, {float(scr.n1.A)}_rt'
            nuc2_info = f'{float(scr.n2.Z)}_rt, {float(scr.n2.A)}_rt'

            if not (scr.n1.dummy or scr.n2.dummy):
                # Scope the screening calculation to avoid multiple definitions of scn_fac.
                of.write(f'\n{self.indent*n_indent}' + '{');

                of.write(f'\n{self.indent*(n_indent+1)}constexpr auto scn_fac = scrn::calculate_screen_factor({nuc1_info}, {nuc2_info});\n\n')

                # Insert a static assert (which will always pass) to require the
                # compiler to evaluate the screen factor at compile time.
                of.write(f'\n{self.indent*(n_indent+1)}static_assert(scn_fac.z1 == {float(scr.n1.Z)}_rt);\n\n');

                of.write(f'\n{self.indent*(n_indent+1)}actual_screen5(pstate, scn_fac, scor, dscor_dt);\n');

                of.write(f'{self.indent*n_indent}' + '}\n\n');

            if scr.name == "he4_he4_he4":
                # we don't need to do anything here, but we want to avoid immediately applying the screening
                pass

            elif scr.name == "he4_he4_he4_dummy":
                # handle the second part of the screening for 3-alpha
                of.write(f'\n{self.indent*n_indent}' + '{');

                of.write(f'\n{self.indent*(n_indent+1)}constexpr auto scn_fac2 = scrn::calculate_screen_factor({nuc1_info}, {nuc2_info});\n\n')

                of.write(f'\n{self.indent*(n_indent+1)}static_assert(scn_fac2.z1 == {float(scr.n1.Z)}_rt);\n\n');

                of.write(f'\n{self.indent*(n_indent+1)}actual_screen5(pstate, scn_fac2, scor2, dscor2_dt);\n');

                of.write(f'\n{self.indent*n_indent}' + '}\n\n');

                # there might be both the forward and reverse 3-alpha
                # if we are doing symmetric screening

                for rr in scr.rates:
                    of.write(f'\n')
                    of.write(f'{self.indent*n_indent}ratraw = rate_eval.screened_rates(k_{rr.fname});\n')
                    of.write(f'{self.indent*n_indent}dratraw_dT = rate_eval.dscreened_rates_dT(k_{rr.fname});\n')
                    of.write(f'{self.indent*n_indent}rate_eval.screened_rates(k_{rr.fname}) *= scor * scor2;\n')
                    of.write(f'{self.indent*n_indent}rate_eval.dscreened_rates_dT(k_{rr.fname}) = ratraw * (scor * dscor2_dt + dscor_dt * scor2) + dratraw_dT * scor * scor2;\n')

            else:
                # there might be several rates that have the same
                # reactants and therefore the same screening applies
                # -- handle them all now

                for rr in scr.rates:
                    of.write(f'\n')
                    of.write(f'{self.indent*n_indent}ratraw = rate_eval.screened_rates(k_{rr.fname});\n')
                    of.write(f'{self.indent*n_indent}dratraw_dT = rate_eval.dscreened_rates_dT(k_{rr.fname});\n')
                    of.write(f'{self.indent*n_indent}rate_eval.screened_rates(k_{rr.fname}) *= scor;\n')
                    of.write(f'{self.indent*n_indent}rate_eval.dscreened_rates_dT(k_{rr.fname}) = ratraw * dscor_dt + dratraw_dT * scor;\n')

            of.write('\n')

        # the C++ screen.H code requires that there be at least 1 screening
        # factor because it statically allocates some arrays, so if we turned
        # off screening, just set num_screen_calls = 1 here.

        self.num_screen_calls = max(1, len(screening_map))


    def _nrat_reaclib(self, n_indent, of):
        # Writes the number of Reaclib rates
        of.write(f'{self.indent*n_indent}const int NrateReaclib = {len(self.reaclib_rates)};\n')

        nreaclib_sets = 0
        for nr in self.reaclib_rates:
            r = self.rates[nr]
            nreaclib_sets = nreaclib_sets + len(r.sets)

        of.write(f'{self.indent*n_indent}const int NumReaclibSets = {nreaclib_sets};\n')

    def _nrat_tabular(self, n_indent, of):
        # Writes the number of tabular rates
        of.write(f'{self.indent*n_indent}const int NrateTabular = {len(self.tabular_rates)};\n')

    def _nrxn(self, n_indent, of):
        for i,r in enumerate(self.rates):
            of.write(f'{self.indent*n_indent}k_{r.fname} = {i+1},\n')
        of.write(f'{self.indent*n_indent}NumRates = k_{self.rates[-1].fname}\n')

    def _ebind(self, n_indent, of):
        for nuc in self.unique_nuclei:
            of.write(f'{self.indent*n_indent}ebind_per_nucleon({nuc.c()}) = {nuc.nucbind}_rt;\n')

    def _screen_add(self, n_indent, of):
        screening_map = self.get_screening_map()
        for scr in screening_map:
            of.write(f'{self.indent*n_indent}add_screening_factor(jscr++, ')
            if not scr.n1.dummy:
                of.write(f'zion[{scr.n1.c()}-1], aion[{scr.n1.c()}-1], ')
            else:
                of.write(f'{float(scr.n1.Z)}_rt, {float(scr.n1.A)}_rt, ')
            if not scr.n2.dummy:
                of.write(f'zion[{scr.n2.c()}-1], aion[{scr.n2.c()}-1]);\n\n')
            else:
                of.write(f'{float(scr.n2.Z)}_rt, {float(scr.n2.A)}_rt);\n\n')

    def _write_reaclib_metadata(self, n_indent, of):
        jset = 0
        for nr in self.reaclib_rates:
            r = self.rates[nr]
            for s in r.sets:
                jset = jset + 1
                for an in s.a:
                    of.write(f'{an}\n')
        j = 1
        for i, r in enumerate(self.rates):
            if i in self.reaclib_rates:
                of.write(f'{j}\n')
                j = j + len(r.sets)

        for i, r in enumerate(self.rates):
            if i in self.reaclib_rates:
                j = len(r.sets)-1
                of.write(f'{j}\n')

    def _table_num(self, n_indent, of):
        of.write(f'{self.indent*n_indent}const int num_tables = {len(self.tabular_rates)};\n')

    def _declare_tables(self, n_indent, of):
        for irate in self.tabular_rates:
            r = self.rates[irate]
            idnt = self.indent*n_indent

            of.write(f'{idnt}extern AMREX_GPU_MANAGED table_t {r.table_index_name}_meta;\n')
            of.write(f'{idnt}extern AMREX_GPU_MANAGED Array3D<Real, 1, {r.table_temp_lines}, 1, {r.table_rhoy_lines}, 1, {r.table_num_vars}> {r.table_index_name}_data;\n')
            of.write(f'{idnt}extern AMREX_GPU_MANAGED Array1D<Real, 1, {r.table_rhoy_lines}> {r.table_index_name}_rhoy;\n')
            of.write(f'{idnt}extern AMREX_GPU_MANAGED Array1D<Real, 1, {r.table_temp_lines}> {r.table_index_name}_temp;\n')
            of.write('\n')

    def _table_declare_meta(self, n_indent, of):
        for irate in self.tabular_rates:
            r = self.rates[irate]
            idnt = self.indent*n_indent

            of.write(f"{idnt}AMREX_GPU_MANAGED table_t {r.table_index_name}_meta;\n")

            of.write(f'{idnt}AMREX_GPU_MANAGED Array3D<Real, 1, {r.table_temp_lines}, 1, {r.table_rhoy_lines}, 1, {r.table_num_vars}> {r.table_index_name}_data;\n')

            of.write(f'{idnt}AMREX_GPU_MANAGED Array1D<Real, 1, {r.table_rhoy_lines}> {r.table_index_name}_rhoy;\n')
            of.write(f'{idnt}AMREX_GPU_MANAGED Array1D<Real, 1, {r.table_temp_lines}> {r.table_index_name}_temp;\n\n')

    def _table_init_meta(self, n_indent, of):
        for irate in self.tabular_rates:
            r = self.rates[irate]
            idnt = self.indent*n_indent
            of.write(f'{idnt}{r.table_index_name}_meta.ntemp = {r.table_temp_lines};\n')
            of.write(f'{idnt}{r.table_index_name}_meta.nrhoy = {r.table_rhoy_lines};\n')
            of.write(f'{idnt}{r.table_index_name}_meta.nvars = {r.table_num_vars};\n')
            of.write(f'{idnt}{r.table_index_name}_meta.nheader = {r.table_header_lines};\n')
            of.write(f'{idnt}{r.table_index_name}_meta.file = "{r.table_file}";\n\n')


            of.write(f'{idnt}init_tab_info({r.table_index_name}_meta, {r.table_index_name}_rhoy, {r.table_index_name}_temp, {r.table_index_name}_data);\n\n')

            of.write('\n')

    def _table_term_meta(self, n_indent, of):
        for irate in self.tabular_rates:
            r = self.rates[irate]

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

    def _table_rates_indices(self, n_indent, of):
        for n,irate in enumerate(self.tabular_rates):
            r = self.rates[irate]
            of.write(f'{self.indent*n_indent}{r.table_index_name}')
            if n != len(self.tabular_rates)-1:
                of.write(', &')
            of.write('\n')

    def _compute_tabular_rates(self, n_indent, of):
        if len(self.tabular_rates) > 0:

            idnt = self.indent*n_indent

            for irate in self.tabular_rates:
                r = self.rates[irate]

                of.write(f'{idnt}tabular_evaluate({r.table_index_name}_meta, {r.table_index_name}_rhoy, {r.table_index_name}_temp, {r.table_index_name}_data,\n')
                of.write(f'{idnt}                 rhoy, state.T, rate, drate_dt, edot_nu);\n')

                of.write(f'{idnt}rate_eval.screened_rates(k_{r.fname}) = rate;\n')
                of.write(f'{idnt}rate_eval.dscreened_rates_dT(k_{r.fname}) = drate_dt;\n')
                of.write(f'{idnt}rate_eval.add_energy_rate(k_{r.fname}) = edot_nu;\n')
                of.write('\n')

    def _ydot(self, n_indent, of):
        # Write YDOT
        for n in self.unique_nuclei:
            of.write(f"{self.indent*n_indent}{self.symbol_rates.name_ydot_nuc}({n.c()}) =\n")
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

        for nr, r in enumerate(self.rates):
            if nr in self.tabular_rates:
                if len(r.reactants) != 1:
                    sys.exit('ERROR: Unknown energy rate corrections for a reaction where the number of reactants is not 1.')
                else:
                    reactant = r.reactants[0]
                    of.write(f'{idnt}enuc += C::Legacy::n_A * {self.symbol_rates.name_y}({reactant.c()}) * rate_eval.add_energy_rate(k_{r.fname});\n')

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
                    of.write(f"{self.indent*n_indent}jac.set({nj.c()}, {ni.c()}, scratch);\n\n")

    def _initial_mass_fractions(self, n_indent, of):
        for i, _ in enumerate(self.unique_nuclei):
            if i == 0:
                of.write(f"{self.indent*n_indent}unit_test.X{i+1} = 1.0\n")
            else:
                of.write(f"{self.indent*n_indent}unit_test.X{i+1} = 0.0\n")

    def _pynucastro_home(self, n_indent, of):
        of.write('{}PYNUCASTRO_HOME := {}\n'.format(self.indent*n_indent,
                                                    os.path.dirname(self.pynucastro_dir)))

    def _secret_code_write(self, n_indent, of):
        of.write(f"{self.indent*n_indent}{self.secret_code}\n")

    def _secret_code_write_reference(self, n_indent, of):
        of.write(f"{self.indent*n_indent}const std::string secret_code_reference = \"{self.secret_code}\";\n")
