"""Support for a pure Fortran reaction network.  These functions will
write the Fortran code necessary to integrate a reaction network
comprised of the rates that are passed in.

"""


import os
import shutil
import re
from collections import OrderedDict
from abc import ABC, abstractmethod

import sympy

from pynucastro.rates import Nucleus
from pynucastro.networks import RateCollection
from pynucastro.nucdata import BindingTable


class BaseFortranNetwork(ABC, RateCollection):
    """Interpret the collection of rates and nuclei and produce the
    Fortran code needed to integrate the network.

    """

    def __init__(self, *args, **kwargs):
        """Initialize the Fortran network.  We take a single argument: a list
        of rate files that will make up the network

        """

        super().__init__(*args, **kwargs)

        # Get the template files for writing this network code
        self.template_files = self._get_template_files()

        # a dictionary of functions to call to handle specific parts
        # of the Fortran template
        self.ftags = OrderedDict()
        self.ftags['<nrates>'] = self._nrates
        self.ftags['<nrat_reaclib>'] = self._nrat_reaclib
        self.ftags['<nrat_tabular>'] = self._nrat_tabular
        self.ftags['<nspec>'] = self._nspec
        self.ftags['<nspec_evolve>'] = self._nspec_evolve
        self.ftags['<network_name>'] = self._network_name
        self.ftags['<nrxn>'] = self._nrxn
        self.ftags['<jion>'] = self._jion
        self.ftags['<spec_names>'] = self._spec_names
        self.ftags['<short_spec_names>'] = self._short_spec_names
        self.ftags['<ebind>'] = self._ebind
        self.ftags['<aion>'] = self._aion
        self.ftags['<aion_inv>'] = self._aion_inv
        self.ftags['<zion>'] = self._zion
        self.ftags['<nion>'] = self._nion
        self.ftags['<screen_add>'] = self._screen_add
        self.ftags['<compute_screening_factors>'] = self._compute_screening_factors
        self.ftags['<write_reaclib_metadata>'] = self._write_reaclib_metadata
        self.ftags['<table_num>'] = self._table_num
        self.ftags['<public_table_indices>'] = self._public_table_indices
        self.ftags['<table_indices>'] = self._table_indices
        self.ftags['<declare_tables>'] = self._declare_tables
        self.ftags['<declare_managed_tables>'] = self._declare_managed_tables
        self.ftags['<table_init_meta>'] = self._table_init_meta
        self.ftags['<table_term_meta>'] = self._table_term_meta
        self.ftags['<table_rates_indices>'] = self._table_rates_indices
        self.ftags['<compute_tabular_rates>'] = self._compute_tabular_rates
        self.ftags['<ydot_declare_scratch>'] = self._ydot_declare_scratch
        self.ftags['<ydot_scratch>'] = self._ydot_scratch
        self.ftags['<ydot>'] = self._ydot
        self.ftags['<enuc_add_energy_rate>'] = self._enuc_add_energy_rate
        self.ftags['<jacnuc_declare_scratch>'] = self._jacnuc_declare_scratch
        self.ftags['<jacnuc_scratch>'] = self._jacnuc_scratch
        self.ftags['<jacnuc>'] = self._jacnuc
        self.ftags['<yinit_nuc>'] = self._yinit_nuc
        self.ftags['<initial_mass_fractions>'] = self._initial_mass_fractions
        self.ftags['<probin_mass_fractions>'] = self._probin_mass_fractions
        self.ftags['<parameters_mass_fractions>'] = self._parameters_mass_fractions
        self.ftags['<final_net_print>'] = self._final_net_print
        self.ftags['<headerline>'] = self._headerline
        self.ftags['<pynucastro_home>'] = self._pynucastro_home
        self.indent = '  '

        self.use_cse = False

        self.float_explicit_num_digits = 17

        self.ydot_out_scratch = None
        self.ydot_out_result  = None
        self.solved_ydot      = False
        self.jac_out_scratch  = None
        self.jac_out_result   = None
        self.jac_null_entries = None
        self.solved_jacobian  = False
        self.symbol_ludict = OrderedDict() # Symbol lookup dictionary

        self.num_screen_calls = None

        # Define these for the particular network
        self.name_rate_data = 'screened_rates'
        self.name_y         = 'Y'
        self.name_ydot      = 'ydot'
        self.name_ydot_nuc  = 'ydot_nuc'
        self.name_jacobian  = 'jac'
        self.name_jacobian_nuc  = 'jac'
        self.name_density   = 'state % rho'
        self.name_electron_fraction = 'state % y_e'
        self.symbol_ludict['__dens__'] = self.name_density
        self.symbol_ludict['__y_e__'] = self.name_electron_fraction

    @abstractmethod
    def _get_template_files(self):
        # This method should be overridden by derived classes
        # to support specific output templates.
        # This method returns a list of strings that are file paths to template files.
        return []

    def ydot_string(self, rate):
        """
        return a string containing the term in a dY/dt equation
        in a reaction network corresponding to this rate for Fortran 90.
        """

        # composition dependence
        Y_string = ""
        for n, r in enumerate(sorted(set(rate.reactants))):
            c = rate.reactants.count(r)
            if c > 1:
                Y_string += self.name_y + f"(j{r})**{c}"
            else:
                Y_string += self.name_y + f"(j{r})"

            if n < len(set(rate.reactants))-1:
                Y_string += " * "

        # density dependence
        if rate.dens_exp == 0:
            dens_string = ""
        elif rate.dens_exp == 1:
            dens_string = f"{self.name_density} * "
        else:
            dens_string = f"{self.name_density}**{rate.dens_exp} * "

        # prefactor
        if not rate.prefactor == 1.0:
            prefactor_string = f"{rate.prefactor:1.14e} * ".replace('e','d')
        else:
            prefactor_string = ""

        return "{}{}{} * {}(i_scor, k_{}) * {}(i_rate, k_{})".format(
            prefactor_string,
            dens_string,
            Y_string,
            self.name_rate_data,
            rate.fname,
            self.name_rate_data,
            rate.fname)

    def ydot_term_symbol(self, rate, y_i):
        """
        return a sympy expression containing this rate's contribution to
        the ydot term for nuclide y_i.
        """
        srate = self.specific_rate_symbol(rate)

        # Check if y_i is a reactant or product
        c_reac = rate.reactants.count(y_i)
        c_prod = rate.products.count(y_i)
        if c_reac == 0 and c_prod == 0:
            # The rate doesn't contribute to the ydot for this y_i
            ydot_sym = float(sympy.sympify(0.0))
        else:
            # y_i appears as a product or reactant
            ydot_sym = (c_prod - c_reac) * srate
        return ydot_sym.evalf(n=self.float_explicit_num_digits)

    def specific_rate_symbol(self, rate):
        """
        return a sympy expression containing the term in a dY/dt equation
        in a reaction network corresponding to this rate.

        Also enter the symbol and substitution in the lookup table.
        """

        # composition dependence
        Y_sym = 1
        for r in sorted(set(rate.reactants)):
            c = rate.reactants.count(r)
            sym_final = self.name_y + f'(j{r})'
            sym_temp  = f'Y__j{r}__'
            self.symbol_ludict[sym_temp] = sym_final
            Y_sym = Y_sym * sympy.symbols(sym_temp)**c

        # density dependence
        dens_sym = sympy.symbols('__dens__')**rate.dens_exp

        # electron fraction if electron capture reaction
        if (rate.weak_type == 'electron_capture' and not rate.tabular):
            y_e_sym = sympy.symbols('__y_e__')
        else:
            y_e_sym = sympy.sympify(1)

        # prefactor
        prefactor_sym = sympy.sympify(1)/sympy.sympify(rate.inv_prefactor)

        # screened rate
        sym_final = self.name_rate_data + f'(k_{rate.fname})'
        sym_temp  = f'NRD__k_{rate.fname}__'
        self.symbol_ludict[sym_temp] = sym_final
        screened_rate_sym = sympy.symbols(sym_temp)

        srate_sym = prefactor_sym * dens_sym * y_e_sym * Y_sym * screened_rate_sym
        return srate_sym

    def fortranify(self, s):
        """
        Given string s, will replace the symbols appearing as keys in
        self.symbol_ludict with their corresponding entries.
        """
        for k in self.symbol_ludict:
            v = self.symbol_ludict[k]
            s = s.replace(k,v)
        if s == '0':
            s = '0.0e0_rt'

        ## Replace all double precision literals with custom real type literals
        # constant type specifier
        const_spec = "_rt"

        # we want to replace any "d" scientific notation with the new style
        # this matches stuff like -1.25d-10, and gives us separate groups for the
        # prefix and exponent.  The [^\w] makes sure a letter isn't right in front
        # of the match (like 'k3d-1'). Alternately, we allow for a match at the start of the string.
        d_re = re.compile(r"([^\w\+\-]|\A)([\+\-0-9.][0-9.]+)[dD]([\+\-]?[0-9]+)", re.IGNORECASE|re.DOTALL)

        # update "d" scientific notation -- allow for multiple constants in a single string
        for dd in d_re.finditer(s):
            prefix = dd.group(2)
            exponent = dd.group(3)
            new_num = f"{prefix}e{exponent}{const_spec}"
            old_num = dd.group(0).strip()
            s = s.replace(old_num, new_num)

        return s

    def jacobian_string(self, rate, ydot_j, y_i):
        """
        return a string containing the term in a jacobian matrix
        in a reaction network corresponding to this rate

        Returns the derivative of the j-th YDOT wrt. the i-th Y
        If the derivative is zero, returns the empty string ''

        ydot_j and y_i are objects of the class 'Nucleus'
        """
        if (ydot_j not in rate.reactants and ydot_j not in rate.products) or \
            y_i not in rate.reactants:
            return ''

        # composition dependence
        Y_string = ""
        for n, r in enumerate(sorted(set(rate.reactants))):
            c = rate.reactants.count(r)
            if y_i == r:
                if c == 1:
                    continue
                if n>0 and n < len(set(rate.reactants))-1:
                    Y_string += "*"
                if c > 2:
                    Y_string += f"{c}*{self.name_y}(j{r})**{c-1}"
                elif c==2:
                    Y_string += f"2*{self.name_y}(j{r})"
            else:
                if n>0 and n < len(set(rate.reactants))-1:
                    Y_string += "*"
                if c > 1:
                    Y_string += f"{self.name_y}(j{r})**{c}"
                else:
                    Y_string += f"{self.name_y}(j{r})"

        # density dependence
        if rate.dens_exp == 0:
            dens_string = ""
        elif rate.dens_exp == 1:
            dens_string = f"{self.name_density} * "
        else:
            dens_string = f"{self.name_density}**{rate.dens_exp} * "

        # prefactor
        if not rate.prefactor == 1.0:
            prefactor_string = f"{rate.prefactor:1.14e} * ".replace('e','d')
        else:
            prefactor_string = ""

        if Y_string=="" and dens_string=="" and prefactor_string=="":
            rstring = "{}{}{}   {}(i_scor, k_{}) * {}(i_rate, k_{})"
        else:
            rstring = "{}{}{} * {}(i_scor, k_{}) * {}(i_rate, k_{})"
        return rstring.format(prefactor_string, dens_string, Y_string,
                              self.name_rate_data, rate.fname,
                              self.name_rate_data, rate.fname)


    def jacobian_term_symbol(self, rate, ydot_j, y_i):
        """
        return a sympy expression containing the term in a jacobian matrix
        in a reaction network corresponding to this rate

        Returns the derivative of the j-th YDOT wrt. the i-th Y
        If the derivative is zero, returns 0.

        ydot_j and y_i are objects of the class 'Nucleus'
        """
        ydot_sym = self.ydot_term_symbol(rate, ydot_j)
        deriv_sym = sympy.symbols(f'Y__j{y_i}__')
        jac_sym = sympy.diff(ydot_sym, deriv_sym)
        symbol_is_null = False
        if jac_sym.equals(0):
            symbol_is_null = True
        return (jac_sym.evalf(n=self.float_explicit_num_digits), symbol_is_null)


    def compose_ydot(self):
        """create the expressions for dYdt for the nuclei, where Y is the
        molar fraction.

        """

        ydot = []
        for n in self.unique_nuclei:
            ydot_sym = float(sympy.sympify(0.0))
            for r in self.nuclei_consumed[n]:
                ydot_sym = ydot_sym + self.ydot_term_symbol(r, n)
            for r in self.nuclei_produced[n]:
                ydot_sym = ydot_sym + self.ydot_term_symbol(r, n)
            ydot.append(ydot_sym)

        if self.use_cse:
            scratch_sym = sympy.utilities.numbered_symbols('scratch_')
            scratch, result = sympy.cse(ydot, symbols=scratch_sym, order='none')

            result_out = []
            for r in result:
                result_out.append(r.evalf(n=self.float_explicit_num_digits))
            scratch_out = []
            for s in scratch:
                scratch_out.append([s[0], s[1].evalf(n=self.float_explicit_num_digits)])
            self.ydot_out_scratch = scratch_out
            self.ydot_out_result  = result_out
        else:
            self.ydot_out_scratch = None
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
                    rsym_add, rsym_add_null = self.jacobian_term_symbol(r, nj, ni)
                    rsym = rsym + rsym_add
                    rsym_is_null = rsym_is_null and rsym_add_null
                for r in self.nuclei_produced[nj]:
                    rsym_add, rsym_add_null = self.jacobian_term_symbol(r, nj, ni)
                    rsym = rsym + rsym_add
                    rsym_is_null = rsym_is_null and rsym_add_null
                jac_sym.append(rsym)
                jac_null.append(rsym_is_null)

        if self.use_cse:
            scratch_sym = sympy.utilities.numbered_symbols('scratch_')
            scratch, result = sympy.cse(jac_sym, symbols=scratch_sym, order='none')

            result_out = []
            for r in result:
                result_out.append(r.evalf(n=self.float_explicit_num_digits))
            scratch_out = []
            for s in scratch:
                scratch_out.append([s[0], s[1].evalf(n=self.float_explicit_num_digits)])
            self.jac_out_scratch = scratch_out
            self.jac_out_result  = result_out
        else:
            self.jac_out_scratch = None
            self.jac_out_result  = jac_sym
        self.jac_null_entries = jac_null
        self.solved_jacobian = True

    def io_open(self, infile, outfile):
        """open the input and output files"""
        try:
            of = open(outfile, "w")
        except:
            raise
        try:
            ifile = open(infile)
        except:
            raise
        return ifile, of

    def io_close(self, infile, outfile):
        """close the input and output files"""
        infile.close()
        outfile.close()

    def fmt_to_rt_f90(self, i):
        """convert a number to custom real type format"""
        return f'{float(i):1.14e}_rt'

    def get_indent_amt(self, l, k):
        """determine the amount of spaces to indent a line"""
        rem = re.match(r'\A'+k+r'\(([0-9]*)\)\Z',l)
        return int(rem.group(1))

    def _write_network(self, use_cse=False):
        """
        This writes the RHS, jacobian and ancillary files for the system of ODEs that
        this network describes, using the template files.
        """
        self.use_cse = use_cse

        # Prepare RHS terms
        if not self.solved_ydot:
            self.compose_ydot()
        if not self.solved_jacobian:
            self.compose_jacobian()

        # Process template files
        for tfile in self.template_files:
            tfile_basename = os.path.basename(tfile)
            outfile    = tfile_basename.replace('.template', '')
            ifile, of = self.io_open(tfile, outfile)
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
            self.io_close(ifile, of)

        # Copy any tables in the network to the current directory
        # if the table file cannot be found, print a warning and continue.
        for i_tab in self.tabular_rates:
            tr = self.rates[i_tab]
            tdir = os.path.dirname(tr.rfile_path)
            if tdir != os.getcwd():
                tdat_file = os.path.join(tdir, tr.table_file)
                if os.path.isfile(tdat_file):
                    shutil.copy(tdat_file, os.getcwd())
                    #Rename *.dat to *.network.dat
                    copy_tdat = os.path.join(os.getcwd(), os.path.basename(tdat_file))
                    renamed_tdat = copy_tdat.replace('.dat', '.network.dat')
                    os.rename(copy_tdat, renamed_tdat)
                else:
                    print(f'WARNING: Table data file {tr.table_file} not found.')

    def _nrates(self, n_indent, of):
        of.write('{}integer, parameter :: nrates = {}\n'.format(
            self.indent*n_indent,
            len(self.rates)))

    def get_screening_map(self):
        """a screening map is just a list of tuples containing the information
        about nuclei pairs for screening:
        (descriptive name of nuclei, nucleus 1, nucleus 2, rate, 1-based index of rate)
        """
        screening_map = []
        for k, r in enumerate(self.rates):
            if r.ion_screen:
                nucs = "_".join([str(q) for q in r.ion_screen])
                in_map = False
                for h, _, _, mrates, krates in screening_map:
                    if h == nucs:
                        # if we already have the reactants, then we
                        # will already be doing the screening factors,
                        # so just append this new rate to the list we
                        # are keeping of the rates where this
                        # screening is needed
                        in_map = True
                        mrates.append(r)
                        krates.append(k+1)
                        break
                if not in_map:
                    # we handle 3-alpha specially -- we actually need 2 screening factors for it
                    if nucs == "he4_he4_he4":
                        # he4 + he4
                        screening_map.append((nucs, r.ion_screen[0], r.ion_screen[1],
                                              [r], [k+1]))
                        # he4 + be8
                        be8 = Nucleus("Be8", dummy=True)
                        screening_map.append((nucs+"_dummy", r.ion_screen[2], be8,
                                              [r], [k+1]))

                    else:
                        screening_map.append((nucs, r.ion_screen[0], r.ion_screen[1],
                                              [r], [k+1]))
        return screening_map

    def _compute_screening_factors(self, n_indent, of):
        screening_map = self.get_screening_map()
        for i, (h, _, _, mrates, krates) in enumerate(screening_map):
            if h == "he4_he4_he4":
                # handle both parts of the 3-alpha screening here
                of.write(f'\n{self.indent*n_indent}call screen5(pstate, {i+1}, scor, dscor_dt, dscor_dd)\n')
                of.write(f'\n{self.indent*n_indent}call screen5(pstate, {i+2}, scor2, dscor2_dt, dscor2_dd)\n')
                of.write(f'{self.indent*n_indent}rate_eval % unscreened_rates(i_scor,{krates[0]}) = scor * scor2\n')
                of.write(f'{self.indent*n_indent}rate_eval % unscreened_rates(i_dscor_dt,{krates[0]}) = scor * dscor2_dt + dscor_dt * scor2\n')

            elif h == "he4_he4_he4_dummy":
                continue

            else:
                of.write(f'\n{self.indent*n_indent}call screen5(pstate, {i+1}, scor, dscor_dt, dscor_dd)\n')
                for _, k in zip(mrates, krates):
                    of.write(f'{self.indent*n_indent}rate_eval % unscreened_rates(i_scor,{k}) = scor\n')
                    of.write(f'{self.indent*n_indent}rate_eval % unscreened_rates(i_dscor_dt,{k}) = dscor_dt\n')

            of.write('\n')

        self.num_screen_calls = len(screening_map)

    def _nrat_reaclib(self, n_indent, of):
        # Writes the number of Reaclib rates
        of.write('{}integer, parameter :: nrat_reaclib = {}\n'.format(
            self.indent*n_indent,
            len(self.reaclib_rates)))

        nreaclib_sets = 0
        for nr in self.reaclib_rates:
            r = self.rates[nr]
            nreaclib_sets = nreaclib_sets + len(r.sets)

        of.write('{}integer, parameter :: number_reaclib_sets = {}\n'.format(
            self.indent*n_indent,
            nreaclib_sets))

    def _nrat_tabular(self, n_indent, of):
        # Writes the number of tabular rates
        of.write('{}integer, parameter :: nrat_tabular = {}\n'.format(
            self.indent*n_indent,
            len(self.tabular_rates)))

    def _nspec(self, n_indent, of):
        of.write('{}integer, parameter :: nspec = {}\n'.format(
            self.indent*n_indent,
            len(self.unique_nuclei)))

    def _nspec_evolve(self, n_indent, of):
        # Evolve all the nuclei at the moment
        of.write('{}integer, parameter :: nspec_evolve = {}\n'.format(
            self.indent*n_indent,
            len(self.unique_nuclei)))

    def _network_name(self, n_indent, of):
        # the name of the network
        of.write('{}character (len=32), parameter :: network_name = "{}"\n'.format(
            self.indent*n_indent,
            "pynucastro"))

    def _jion(self, n_indent, of):
        for i,nuc in enumerate(self.unique_nuclei):
            of.write('{}integer, parameter :: j{}   = {}\n'.format(
                self.indent*n_indent, nuc, i+1))

    def _spec_names(self, n_indent, of):
        for nuc in self.unique_nuclei:
            of.write('{}spec_names(j{})   = "{}"\n'.format(
                self.indent*n_indent, nuc, nuc.spec_name))

    def _short_spec_names(self, n_indent, of):
        for nuc in self.unique_nuclei:
            of.write('{}short_spec_names(j{})   = "{}"\n'.format(
                self.indent*n_indent, nuc, nuc.short_spec_name))

    def _nrxn(self, n_indent, of):
        for i,r in enumerate(self.rates):
            of.write('{}integer, parameter :: k_{}   = {}\n'.format(
                self.indent*n_indent, r.fname, i+1))

    def _ebind(self, n_indent, of):
        bintable = BindingTable()
        for nuc in self.unique_nuclei:
            nuc_in_table = bintable.get_nuclide(n=nuc.N, z=nuc.Z)
            str_nucbind = self.fmt_to_rt_f90(nuc_in_table.nucbind)
            of.write('{}ebind_per_nucleon(j{})   = {}\n'.format(
                self.indent*n_indent, nuc, str_nucbind))

    def _aion(self, n_indent, of):
        for nuc in self.unique_nuclei:
            of.write('{}aion(j{})   = {}\n'.format(
                self.indent*n_indent,
                nuc,
                self.fmt_to_rt_f90(nuc.A)))

    def _aion_inv(self, n_indent, of):
        for nuc in self.unique_nuclei:
            of.write('{}aion_inv(j{})   = 1.0_rt/{}\n'.format(
                self.indent*n_indent,
                nuc,
                self.fmt_to_rt_f90(nuc.A)))

    def _zion(self, n_indent, of):
        for nuc in self.unique_nuclei:
            of.write('{}zion(j{})   = {}\n'.format(
                self.indent*n_indent,
                nuc,
                self.fmt_to_rt_f90(nuc.Z)))

    def _nion(self, n_indent, of):
        for nuc in self.unique_nuclei:
            of.write('{}nion(j{})   = {}\n'.format(
                self.indent*n_indent,
                nuc,
                self.fmt_to_rt_f90(nuc.N)))

    def _screen_add(self, n_indent, of):
        screening_map = self.get_screening_map()
        for _, n1, n2, _, _ in screening_map:
            of.write(f'{self.indent*n_indent}call add_screening_factor(')
            if not n1.dummy:
                of.write(f'zion(j{n1}), aion(j{n1}), &\n')
            else:
                of.write(f'{float(n1.Z)}_rt, {float(n1.A)}_rt), &\n')
            if not n2.dummy:
                of.write(f'{self.indent*(n_indent+1)}zion(j{n2}), aion(j{n2}))\n\n')
            else:
                of.write(f'{self.indent*(n_indent+1)}{float(n2.Z)}_rt, {float(n2.A)}_rt)\n\n')

    def _write_reaclib_metadata(self, n_indent, of):
        jset = 0
        for nr in self.reaclib_rates:
            r = self.rates[nr]
            for _ in r.sets:
                jset = jset + 1
                for an in a:
                    of.write(f'{self.fmt_to_dp_f90(an)}\n')
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
        of.write('{}integer, parameter :: num_tables   = {}\n'.format(
            self.indent*n_indent, len(self.tabular_rates)))

    def _public_table_indices(self, n_indent, of):
        for irate in self.tabular_rates:
            r = self.rates[irate]
            of.write(f'{self.indent*n_indent}public {r.table_index_name}\n')

    def _table_indices(self, n_indent, of):
        for n,irate in enumerate(self.tabular_rates):
            r = self.rates[irate]
            of.write('{}integer, parameter :: {}   = {}\n'.format(
                self.indent*n_indent, r.table_index_name, n+1))

    def _declare_tables(self, n_indent, of):
        for irate in self.tabular_rates:
            r = self.rates[irate]
            of.write('{}real(rt), allocatable :: rate_table_{}(:,:,:), rhoy_table_{}(:), temp_table_{}(:)\n'.format(
                self.indent*n_indent, r.table_index_name, r.table_index_name, r.table_index_name))
            of.write('{}integer, allocatable  :: num_rhoy_{}, num_temp_{}, num_vars_{}\n'.format(
                self.indent*n_indent, r.table_index_name, r.table_index_name, r.table_index_name))
            of.write('{}character(len=50)     :: rate_table_file_{}\n'.format(
                self.indent*n_indent, r.table_index_name))
            of.write('{}integer               :: num_header_{}\n'.format(
                self.indent*n_indent, r.table_index_name))
            of.write('\n')

    def _declare_managed_tables(self, n_indent, of):
        for irate in self.tabular_rates:
            r = self.rates[irate]
            of.write('{}attributes(managed) :: rate_table_{}, rhoy_table_{}, temp_table_{}\n'.format(
                self.indent*n_indent, r.table_index_name, r.table_index_name, r.table_index_name))
            of.write('{}attributes(managed) :: num_rhoy_{}, num_temp_{}, num_vars_{}\n'.format(
                self.indent*n_indent, r.table_index_name, r.table_index_name, r.table_index_name))
            of.write('\n')

    def _table_init_meta(self, n_indent, of):
        for irate in self.tabular_rates:
            r = self.rates[irate]

            of.write('{}allocate(num_temp_{})\n'.format(
                self.indent*n_indent, r.table_index_name))

            of.write('{}allocate(num_rhoy_{})\n'.format(
                self.indent*n_indent, r.table_index_name))

            of.write('{}allocate(num_vars_{})\n'.format(
                self.indent*n_indent, r.table_index_name))

            of.write('{}num_temp_{} = {}\n'.format(
                self.indent*n_indent, r.table_index_name, r.table_temp_lines))

            of.write('{}num_rhoy_{} = {}\n'.format(
                self.indent*n_indent, r.table_index_name, r.table_rhoy_lines))

            of.write('{}num_vars_{} = {}\n'.format(
                self.indent*n_indent, r.table_index_name, r.table_num_vars))

            of.write('{}num_header_{} = {}\n'.format(
                self.indent*n_indent, r.table_index_name, r.table_header_lines))

            of.write('{}rate_table_file_{} = trim("{}")\n'.format(
                self.indent*n_indent, r.table_index_name, r.table_file))

            of.write('{}allocate(rate_table_{}(num_temp_{}, num_rhoy_{}, num_vars_{}))\n'.format(
                self.indent*n_indent, r.table_index_name, r.table_index_name, r.table_index_name, r.table_index_name))

            of.write('{}allocate(rhoy_table_{}(num_rhoy_{}))\n'.format(
                self.indent*n_indent, r.table_index_name, r.table_index_name))

            of.write('{}allocate(temp_table_{}(num_temp_{}))\n'.format(
                self.indent*n_indent, r.table_index_name, r.table_index_name))

            of.write('{}call init_tab_info(rate_table_{}, rhoy_table_{}, temp_table_{}, num_rhoy_{}, num_temp_{}, num_vars_{}, rate_table_file_{}, num_header_{})\n'.format(
                self.indent*n_indent, r.table_index_name, r.table_index_name, r.table_index_name, r.table_index_name,
                r.table_index_name, r.table_index_name, r.table_index_name, r.table_index_name))

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
            of.write(f'{self.indent*n_indent}! Calculate tabular rates\n')
            for n, irate in enumerate(self.tabular_rates):
                r = self.rates[irate]
                of.write(f'{self.indent*n_indent}call tabular_evaluate(rate_table_{r.table_index_name}, rhoy_table_{r.table_index_name}, temp_table_{r.table_index_name}, &\n')
                of.write(f'{self.indent*n_indent}                      num_rhoy_{r.table_index_name}, num_temp_{r.table_index_name}, num_vars_{r.table_index_name}, &\n')
                of.write(f'{self.indent*n_indent}                      rhoy, state % T, rate, drate_dt, edot_nu)\n')
                of.write(f'{self.indent*n_indent}rate_eval % unscreened_rates(i_rate,{n+1+len(self.reaclib_rates)}) = rate\n')
                of.write(f'{self.indent*n_indent}rate_eval % unscreened_rates(i_drate_dt,{n+1+len(self.reaclib_rates)}) = drate_dt\n')
                of.write(f'{self.indent*n_indent}rate_eval % add_energy_rate({n+1})  = edot_nu\n')
                of.write('\n')

    def _ydot_declare_scratch(self, n_indent, of):
        # Declare scratch variables
        if self.use_cse:
            for si in self.ydot_out_scratch:
                siname = si[0]
                of.write(f'{self.indent*n_indent}double precision :: {siname}\n')

    def _ydot_scratch(self, n_indent, of):
        # Assign scratch variables
        if self.use_cse:
            for si in self.ydot_out_scratch:
                siname = si[0]
                sivalue = self.fortranify(sympy.fcode(si[1], precision=15,
                                                      source_format='free',
                                                      standard=95))
                of.write(f'{self.indent*n_indent}{siname} = {sivalue}\n')

    def _ydot(self, n_indent, of):
        # Write YDOT
        for i, n in enumerate(self.unique_nuclei):
            sol_value = self.fortranify(sympy.fcode(self.ydot_out_result[i], precision=15,
                                                    source_format='free',
                                                    standard=95))
            of.write('{}{}(j{}) = ( &\n'.format(self.indent*n_indent,
                                                self.name_ydot_nuc, n))
            of.write(f"{self.indent*(n_indent+1)}{sol_value} &\n")
            of.write(f"{self.indent*n_indent}   )\n\n")

    def _enuc_add_energy_rate(self, n_indent, of):
        # Add tabular per-reaction neutrino energy generation rates to the energy generation rate
        # (not thermal neutrinos)
        for nr, r in enumerate(self.rates):
            if nr in self.tabular_rates:
                if len(r.reactants) != 1:
                    print('ERROR: Unknown energy rate corrections for a reaction where the number of reactants is not 1.')
                    exit()
                else:
                    reactant = r.reactants[0]
                    of.write('{}enuc = enuc + N_AVO * {}(j{}) * rate_eval % add_energy_rate({})\n'.format(
                        self.indent*n_indent, self.name_y, reactant, r.table_index_name))

    def _jacnuc_declare_scratch(self, n_indent, of):
        # Declare scratch variables
        if self.use_cse:
            for si in self.jac_out_scratch:
                siname = si[0]
                of.write(f'{self.indent*n_indent}double precision :: {siname}\n')

    def _jacnuc_scratch(self, n_indent, of):
        # Assign scratch variables
        if self.use_cse:
            for si in self.jac_out_scratch:
                siname = si[0]
                sivalue = self.fortranify(sympy.fcode(si[1], precision=15,
                                                      source_format='free',
                                                      standard=95))
                of.write(f'{self.indent*n_indent}{siname} = {sivalue}\n')

    def _jacnuc(self, n_indent, of):
        # now make the Jacobian
        n_unique_nuclei = len(self.unique_nuclei)
        for jnj, nj in enumerate(self.unique_nuclei):
            for ini, ni in enumerate(self.unique_nuclei):
                jac_idx = n_unique_nuclei*jnj + ini
                if not self.jac_null_entries[jac_idx]:
                    jvalue = self.fortranify(sympy.fcode(self.jac_out_result[jac_idx],
                                                         precision=15,
                                                         source_format='free',
                                                         standard=95))
                    of.write(f"{self.indent*(n_indent)}scratch = (&\n")
                    of.write(f"{self.indent*(n_indent+1)}{jvalue} &\n")
                    of.write(f"{self.indent*n_indent}   )\n")
                    of.write("{}call set_jac_entry({}, j{}, j{}, scratch)\n\n".format(
                        self.indent*n_indent, self.name_jacobian, nj, ni))

    def _yinit_nuc(self, n_indent, of):
        for n in self.unique_nuclei:
            of.write(f"{self.indent*n_indent}state_in % xn(j{n}) = initial_mass_fraction_{n}\n")

    def _initial_mass_fractions(self, n_indent, of):
        for n in self.unique_nuclei:
            of.write(f"{self.indent*n_indent}initial_mass_fraction_{n} = 0.0d0\n")

    def _probin_mass_fractions(self, n_indent, of):
        num_unique_nuclei = len(self.unique_nuclei)
        for j, n in enumerate(self.unique_nuclei):
            of.write("{self.indent*n_indent}initial_mass_fraction_{n}")
            if j < num_unique_nuclei - 1:
                of.write(", &\n")

    def _parameters_mass_fractions(self, n_indent, of):
        for n in self.unique_nuclei:
            of.write("{}initial_mass_fraction_{}          real          0.0d0\n".format(
                self.indent*n_indent, n))

    def _final_net_print(self, n_indent, of):
        for n in self.unique_nuclei:
            of.write(f"{self.indent*n_indent}write(*,'(A,ES25.14)') '{n}: ', history % X(j{n}, end_index)\n")

    def _headerline(self, n_indent, of):
        of.write(f'{self.indent*n_indent}write(2, fmt=hfmt) ')
        of.write("'Time', ")
        for nuc in self.unique_nuclei:
            of.write(f"'Y_{nuc}', ")
        of.write("'E_nuc'\n")

    def _pynucastro_home(self, n_indent, of):
        of.write('{}PYNUCASTRO_HOME := {}\n'.format(self.indent*n_indent,
                                                    os.path.dirname(self.pynucastro_dir)))
