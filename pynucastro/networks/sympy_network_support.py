"""This is a module that interprets the rates, ydots, and Jacobian
through sympy"""

import re

import sympy

from pynucastro.rates import TabularRate


class SympyRates:

    def __init__(self):

        self.symbol_ludict = {}  # Symbol lookup dictionary
        self._ydot_term_cache = {}

        self.name_density = 'state.rho'
        self.name_electron_fraction = 'state.y_e'

        # Define these for the particular network
        self.name_rate_data = 'screened_rates'
        self.name_y = 'Y'
        self.name_ydot = 'ydot'
        self.name_ydot_nuc = 'ydot_nuc'
        self.name_jacobian = 'jac'
        self.name_jacobian_nuc = 'jac'
        self.symbol_ludict['__dens__'] = self.name_density
        self.symbol_ludict['__y_e__'] = self.name_electron_fraction

        self.float_explicit_num_digits = 17

    def ydot_term_symbol(self, rate, y_i):
        """
        return a sympy expression containing this rate's contribution to
        the ydot term for nuclide y_i.
        """
        key = (rate.cname(), y_i)
        if key in self._ydot_term_cache:
            return self._ydot_term_cache[key]
        srate = self.specific_rate_symbol(rate)

        # Check if y_i is a reactant or product
        c_reac = rate.reactants.count(y_i)
        c_prod = rate.products.count(y_i)
        if c_reac == 0 and c_prod == 0:
            # The rate doesn't contribute to the ydot for this y_i
            ydot_sym = float(sympy.sympify(0.0))
        else:
            # y_i appears as a product or reactant
            if rate.stoichiometry:
                # custom stoichiometry is not supported when
                # the same nucleus appears on both the left
                # and right side of the reaction
                assert c_reac == 0 or c_prod == 0
                if c_reac > 0:
                    c = -rate.stoichiometry[y_i]
                else:
                    c = rate.stoichiometry[y_i]
            else:
                c = c_prod - c_reac
            ydot_sym = c * srate
        result = ydot_sym.evalf(n=self.float_explicit_num_digits)
        self._ydot_term_cache[key] = result
        return result

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
            sym_final = f'{self.name_y}({r.cindex()})'
            sym_temp = f'Y__j{r}__'
            self.symbol_ludict[sym_temp] = sym_final
            Y_sym = Y_sym * sympy.symbols(sym_temp)**c

        # density dependence
        dens_sym = sympy.symbols('__dens__')**rate.dens_exp

        # electron fraction if electron capture reaction
        y_e_sym = sympy.sympify(1)
        if not isinstance(rate, TabularRate):
            if rate.weak_type == 'electron_capture' and not rate.tabular:
                y_e_sym = sympy.symbols('__y_e__')

        # prefactor
        prefactor_sym = sympy.sympify(1)/sympy.sympify(rate.inv_prefactor)

        # screened rate
        sym_final = self.name_rate_data + f'(k_{rate.cname()})'
        sym_temp = f'NRD__k_{rate.cname()}__'
        self.symbol_ludict[sym_temp] = sym_final
        screened_rate_sym = sympy.symbols(sym_temp)

        srate_sym = prefactor_sym * dens_sym * y_e_sym * Y_sym * screened_rate_sym
        return srate_sym

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

    def cxxify(self, s):
        """
        Given string s, will replace the symbols appearing as keys in
        self.symbol_ludict with their corresponding entries.
        """
        for k, v in self.symbol_ludict.items():
            s = s.replace(k, v)
        if s == '0':
            s = '0.0e0'

        # Replace all double precision literals with custom real type
        # literals
        # constant type specifier
        const_spec = "_rt"

        # we want append any "e" scientific notation with "_rt".  This
        # matches stuff like -1.25d-10, and gives us separate groups
        # for the prefix and exponent.  The [^\w] makes sure a letter
        # isn't right in front of the match (like
        # 'k3d-1'). Alternately, we allow for a match at the start of
        # the string.
        e_re = re.compile(r"([^\w\+\-]|\A)([\+\-0-9.][0-9.]+)[eE]([\+\-]?[0-9]+)", re.IGNORECASE | re.DOTALL)

        # update "d" scientific notation -- allow for multiple
        # constants in a single string
        for ee in e_re.finditer(s):
            old_num = ee.group(0).strip()
            s = s.replace(old_num, f"{old_num}{const_spec}")

        return s
