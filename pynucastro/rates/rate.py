"""
Classes and methods to interface with files storing rate data.
"""

import os
import re
import io
import numpy as np
import matplotlib.pyplot as plt
import numba

try:
    from numba.experimental import jitclass
except ImportError:
    from numba import jitclass

from pynucastro.nucdata import UnidentifiedElement, PeriodicTable, PartitionFunctionCollection, BindingTable, SpinTable

_pynucastro_dir = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
_pynucastro_rates_dir = os.path.join(_pynucastro_dir, 'library')
_pynucastro_tabular_dir = os.path.join(_pynucastro_rates_dir, 'tabular')


#read the spin table once and store it at the module-level
_spin_table = SpinTable(set_double_gs=False)

# read the binding energy table once and store it at the module-level
_binding_table = BindingTable()


_pcollection = PartitionFunctionCollection(use_high_temperatures = True, use_set='frdm')

def _find_rate_file(ratename):
    """locate the Reaclib or tabular rate or library file given its name.  Return
    None if the file cannot be located, otherwise return its path."""

    # check to see if the rate file is in the working dir or
    # is already the full path
    x = ratename
    if os.path.isfile(x):
        return os.path.realpath(x)

    # check to see if the rate file is in pynucastro/library
    x = os.path.join(_pynucastro_rates_dir, ratename)
    if os.path.isfile(x):
        return os.path.realpath(x)

    # check to see if the rate file is in pynucastro/library/tabular
    x = os.path.join(_pynucastro_tabular_dir, ratename)
    if os.path.isfile(x):
        return os.path.realpath(x)

    # notify user we can't find the file
    raise Exception(f'File {ratename} not found in the working directory, {_pynucastro_rates_dir}, or {_pynucastro_tabular_dir}')



Tfactor_spec = [
('T9', numba.float64),
('T9i', numba.float64),
('T913', numba.float64),
('T913i', numba.float64),
('T953', numba.float64),
('lnT9', numba.float64)
]

@jitclass(Tfactor_spec)
class Tfactors:
    """ precompute temperature factors for speed

    :param float T: input temperature (Kelvin)
    :var T9:    T / 1.e9 K
    :var T9i:   1.0 / T9
    :var T913i  1.0 / T9 ** (1/3)
    :var T913   T9 ** (1/3)
    :var T953   T9 ** (5/3)
    :var lnT9   log(T9)
    """

    def __init__(self, T):
        """ return the Tfactors object.  Here, T is temperature in Kelvin """
        self.T9 = T/1.e9
        self.T9i = 1.0/self.T9
        self.T913i = self.T9i**(1./3.)
        self.T913 = self.T9**(1./3.)
        self.T953 = self.T9**(5./3.)
        self.lnT9 = np.log(self.T9)

class SingleSet:
    """ a set in Reaclib is one piece of a rate, in the form

        lambda = exp[ a_0 + sum_{i=1}^5  a_i T_9**(2i-5)/3  + a_6 log T_9]

    A single rate in Reaclib can be composed of multiple sets

    :param a: the coefficients of the exponential fit
    :param labelprops: a collection of flags that classify a ReacLib rate

    """

    def __init__(self, a, labelprops=None):
        """here a is iterable (e.g., list or numpy array), storing the
           coefficients, a0, ..., a6

        """
        self.a = a
        self.labelprops = labelprops
        self.label = None
        self.resonant = None
        self.weak = None
        self.reverse = None

        self._update_label_properties()

    def _update_label_properties(self):
        """ Set label and flags indicating Set is resonant,
            weak, or reverse. """
        assert isinstance(self.labelprops, str)
        try:
            assert len(self.labelprops) == 6
        except:
            raise
        else:
            self.label = self.labelprops[0:4]
            self.resonant = self.labelprops[4] == 'r'
            self.weak = self.labelprops[4] == 'w'
            self.reverse = self.labelprops[5] == 'v'

    def __eq__(self, other):
        """ Determine whether two SingleSet objects are equal to each other. """
        x = True

        for ai, aj in zip(self.a, other.a):
            x = x and (ai == aj)

        x = x and (self.label == other.label)
        x = x and (self.resonant == other.resonant)
        x = x and (self.weak == other.weak)
        x = x and (self.reverse == other.reverse)
        return x

    def f(self):
        """
        return a function for this set -- note: Tf here is a Tfactors
        object
        """
        return lambda tf: np.exp(self.a[0] +
                                 self.a[1]*tf.T9i +
                                 self.a[2]*tf.T913i +
                                 self.a[3]*tf.T913 +
                                 self.a[4]*tf.T9 +
                                 self.a[5]*tf.T953 +
                                 self.a[6]*tf.lnT9)

    def set_string(self, prefix="set", plus_equal=False):
        """
        return a string containing the python code for this set
        """
        if plus_equal:
            string = f"{prefix} += np.exp( "
        else:
            string = f"{prefix} = np.exp( "
        string += f" {self.a[0]}"
        if not self.a[1] == 0.0:
            string += f" + {self.a[1]}*tf.T9i"
        if not self.a[2] == 0.0:
            string += f" + {self.a[2]}*tf.T913i"
        if not self.a[3] == 0.0:
            string += f" + {self.a[3]}*tf.T913"
        if not (self.a[4] == 0.0 and self.a[5] == 0.0 and self.a[6] == 0.0):
            string += "\n{}         ".format(len(prefix)*" ")
        if not self.a[4] == 0.0:
            string += f" + {self.a[4]}*tf.T9"
        if not self.a[5] == 0.0:
            string += f" + {self.a[5]}*tf.T953"
        if not self.a[6] == 0.0:
            string += f" + {self.a[6]}*tf.lnT9"
        string += ")"
        return string


class UnsupportedNucleus(BaseException):
    def __init__(self):
        return


class Nucleus:
    """
    a nucleus that participates in a reaction -- we store it in a
    class to hold its properties, define a sorting, and give it a
    pretty printing string.

    :var Z:               atomic number
    :var N:               neutron number
    :var A:               atomic mass
    :var nucbind:         nuclear binding energy (MeV / nucleon)
    :var short_spec_name: nucleus abbrevation (e.g. "he4")
    :var caps_name:       capitalized short species name (e.g. "He4")
    :var el:              element name (e.g. "he")
    :var pretty:          LaTeX formatted version of the nucleus name

    """
    def __init__(self, name, dummy=False):
        name = name.lower()
        self.raw = name

        # a dummy nucleus is one that we can use where a nucleus is needed
        # but it is not considered to be part of the network
        self.dummy = dummy

        # element symbol and atomic weight
        if name == "p":
            self.el = "h"
            self.A = 1
            self.short_spec_name = "h1"
            self.caps_name = "H1"
        elif name == "d":
            self.el = "h"
            self.A = 2
            self.short_spec_name = "h2"
            self.caps_name = "H2"
        elif name == "t":
            self.el = "h"
            self.A = 3
            self.short_spec_name = "h3"
            self.caps_name = "H3"
        elif name == "a":
            #this is a convenience, enabling the use of a commonly-used alias:
            #    He4 --> \alpha --> "a" , e.g. c12(a,g)o16
            self.el ="he"
            self.A = 4
            self.short_spec_name = "he4"
            self.raw = "he4"
            self.caps_name = "He4"
        elif name == "n":
            self.el = "n"
            self.A = 1
            self.Z = 0
            self.N = 1
            self.short_spec_name = "n"
            self.spec_name = "neutron"
            self.pretty = fr"\mathrm{{{self.el}}}"
            self.caps_name = "N"
        else:
            e = re.match(r"([a-zA-Z]*)(\d*)", name)
            self.el = e.group(1).title()  # chemical symbol
            assert self.el
            try:
                self.A = int(e.group(2))
            except:
                if (name.strip() == 'al-6' or
                    name.strip() == 'al*6'):
                    raise UnsupportedNucleus()
                else:
                    raise
            assert self.A >= 0
            self.short_spec_name = name
            self.caps_name = name.capitalize()

        # set the number of spin states
        try:
            self.spin_states = _spin_table.get_spin_nuclide(self.short_spec_name).spin_states
        except NotImplementedError:
            self.spin_states = None

        # use lowercase element abbreviation regardless the case of the input
        self.el = self.el.lower()

        # set a partition function object to every nucleus
        self.partition_function = _pcollection.get_partition_function(self.short_spec_name)

        # atomic number comes from periodic table
        if name != "n":
            try:
                i = PeriodicTable.lookup_abbreviation(self.el)
            except UnidentifiedElement:
                print(f'Could not identify element: {self.el}')
                raise
            except:
                raise
            else:
                self.Z = i.Z
                assert isinstance(self.Z, int)
                assert self.Z >= 0
                self.N = self.A - self.Z
                assert isinstance(self.N, int)
                assert self.N >= 0

                # long name
                self.spec_name = f'{i.name}-{self.A}'

                # latex formatted style
                self.pretty = fr"{{}}^{{{self.A}}}\mathrm{{{self.el.capitalize()}}}"

        try:
            self.nucbind = _binding_table.get_nuclide(n=self.N, z=self.Z).nucbind
        except NotImplementedError:
            # the binding energy table doesn't know about this nucleus
            self.nucbind = None

    def __repr__(self):
        return self.raw

    def __hash__(self):
        return hash((self.Z, self.A))

    def c(self):
        """return the name capitalized"""
        return self.caps_name

    def __eq__(self, other):
        if isinstance(other, Nucleus):
            return self.el == other.el and \
               self.Z == other.Z and self.A == other.A
        if isinstance(other, tuple):
            return (self.Z, self.A) == other
        return NotImplemented

    def __lt__(self, other):
        if not self.Z == other.Z:
            return self.Z < other.Z
        return self.A < other.A

class Rate:
    """A single reaction rate.  Currently, this can be a
    Reaclib rate, which can be composed of multiple sets, or a tabulated
    electron capture rate."""
    def __init__(self, rfile=None, rfile_path=None, chapter=None, original_source=None,
                 reactants=None, products=None, sets=None, labelprops=None, Q=None):
        """ rfile can be either a string specifying the path to a rate file or
        an io.StringIO object from which to read rate information. """

        self.rfile_path = rfile_path
        self.rfile = None

        if type(rfile) == str:
            self.rfile_path = _find_rate_file(rfile)
            self.rfile = os.path.basename(rfile)

        self.chapter = chapter    # the Reaclib chapter for this reaction
        self.original_source = original_source   # the contents of the original rate file
        self.fname = None

        if reactants:
            self.reactants = reactants
        else:
            self.reactants = []

        if products:
            self.products = products
        else:
            self.products = []

        if sets:
            self.sets = sets
        else:
            self.sets = []

        self.labelprops = labelprops

        self.label = None
        self.resonant = None
        self.resonance_combined = None
        self.weak = None
        self.weak_type = None
        self.reverse = None
        self.tabular = None

        self.Q = Q

        if type(rfile) == str:
            # read in the file, parse the different sets and store them as
            # SingleSet objects in sets[]
            f = open(self.rfile_path)
        elif type(rfile) == io.StringIO:
            # Set f to the io.StringIO object
            f = rfile
        else:
            f = None

        if f:
            self._read_from_file(f)
            f.close()
        else:
            self._set_label_properties()

        self._set_rhs_properties()
        self._set_screening()
        self._set_print_representation()

        if self.tabular:
            self.get_tabular_rate()

    def __repr__(self):
        return self.string

    def __hash__(self):
        return hash(self.__repr__())

    def __eq__(self, other):
        """ Determine whether two Rate objects are equal.
        They are equal if they contain identical reactants and products and
        if they contain the same SingleSet sets and if their chapters are equal."""
        x = True

        x = x and (self.chapter == other.chapter)
        x = x and (self.reactants == other.reactants)
        x = x and (self.products == other.products)
        x = x and (len(self.sets) == len(other.sets))

        for si in self.sets:
            scomp = False
            for sj in other.sets:
                if si == sj:
                    scomp = True
                    break
            x = x and scomp

        return x

    def __lt__(self, other):
        """sort such that lightest reactants come first, and then look at products"""

        # this sort will make two nuclei with the same A be in order of Z
        # (assuming there are no nuclei with A > 999
        # we want to compare based on the heaviest first, so we reverse

        self_react_sorted = sorted(self.reactants, key=lambda x: 1000*x.A + x.Z, reverse=True)
        other_react_sorted = sorted(other.reactants, key=lambda x: 1000*x.A + x.Z, reverse=True)

        if self_react_sorted != other_react_sorted:
            # reactants are different, so now we can check them
            for srn, orn in zip(self_react_sorted, other_react_sorted):
                if not srn == orn:
                    return srn < orn
        else:
            # reactants are the same, so consider products
            self_prod_sorted = sorted(self.products, key=lambda x: 1000*x.A + x.Z, reverse=True)
            other_prod_sorted = sorted(other.products, key=lambda x: 1000*x.A + x.Z, reverse=True)

            for spn, opn in zip(self_prod_sorted, other_prod_sorted):
                if not spn == opn:
                    return spn < opn

        # if we made it here, then the rates are the same
        return True

    def __add__(self, other):
        """Combine the sets of two Rate objects if they describe the same
           reaction. Must be Reaclib rates."""
        assert self.reactants == other.reactants
        assert self.products == other.products
        assert self.chapter == other.chapter
        assert isinstance(self.chapter, int)
        assert self.label == other.label
        assert self.weak == other.weak
        assert self.weak_type == other.weak_type
        assert self.tabular == other.tabular
        assert self.reverse == other.reverse

        if self.resonant != other.resonant:
            self._labelprops_combine_resonance()
        new_rate = Rate(chapter=self.chapter,
                        original_source='\n'.join([self.original_source,
                                                   other.original_source]),
                        reactants=self.reactants,
                        products=self.products,
                        sets=self.sets + other.sets,
                        labelprops=self.labelprops,
                        Q=self.Q)
        return new_rate

    def _set_label_properties(self, labelprops=None):
        """ Calls _update_resonance_combined and then
            _update_label_properties. """
        if labelprops:
            self.labelprops = labelprops

        # Update labelprops based on the Sets in this Rate
        # to set the resonance_combined flag properly
        self._update_resonance_combined()
        self._update_label_properties()

    def _update_resonance_combined(self):
        """ Checks the Sets in this Rate and updates the
            resonance_combined flag as well as
            self.labelprops[4] """
        sres = [s.resonant for s in self.sets]
        if True in sres and False in sres:
            self._labelprops_combine_resonance()
        else:
            self.resonance_combined = False

    def _labelprops_combine_resonance(self):
        """ Update self.labelprops[4] = 'c'.
            Also set the resonance_combined flag. """
        llp = list(self.labelprops)
        llp[4] = 'c'
        self.labelprops = ''.join(llp)
        self.resonance_combined = True

    def _update_label_properties(self):
        """ Set label and flags indicating Rate is resonant,
            weak, or reverse. """
        assert isinstance(self.labelprops, str)
        try:
            assert len(self.labelprops) == 6
        except:
            assert self.labelprops == 'tabular'
            self.label = 'tabular'
            self.resonant = False
            self.resonance_combined = False
            self.weak = False # The tabular rate might or might not be weak
            self.weak_type = None
            self.reverse = False
            self.tabular = True
        else:
            self.label = self.labelprops[0:4]
            self.resonant = self.labelprops[4] == 'r'
            self.weak = self.labelprops[4] == 'w'
            if self.weak:
                if self.label.strip() == 'ec' or self.label.strip() == 'bec':
                    self.weak_type = 'electron_capture'
                else:
                    self.weak_type = self.label.strip().replace('+','_pos_').replace('-','_neg_')
            else:
                self.weak_type = None
            self.reverse = self.labelprops[5] == 'v'
            self.tabular = False

    def _read_from_file(self, f):
        """ given a file object, read rate data from the file. """
        lines = f.readlines()
        f.close()

        self.original_source = "".join(lines)

        # first line is the chapter
        self.chapter = lines[0].strip()
        # catch table prescription
        if self.chapter != "t":
            self.chapter = int(self.chapter)

        # remove any blank lines
        set_lines = [l for l in lines[1:] if not l.strip() == ""]

        if self.chapter == "t":
            # e1 -> e2, Tabulated
            s1 = set_lines.pop(0)
            s2 = set_lines.pop(0)
            s3 = set_lines.pop(0)
            s4 = set_lines.pop(0)
            s5 = set_lines.pop(0)
            f = s1.split()
            try:
                self.reactants.append(Nucleus(f[0]))
                self.products.append(Nucleus(f[1]))
            except:
                print(f'Nucleus objects not be identified in {self.original_source}')
                raise

            self.table_file = s2.strip()
            self.table_header_lines = int(s3.strip())
            self.table_rhoy_lines = int(s4.strip())
            self.table_temp_lines = int(s5.strip())
            self.table_num_vars = 6 # Hard-coded number of variables in tables for now.
            self.table_index_name = f'j_{self.reactants[0]}_{self.products[0]}'
            self.labelprops = 'tabular'
            self._set_label_properties()

        else:
            # the rest is the sets
            first = 1
            while len(set_lines) > 0:
                # check for a new chapter id in case of Reaclib v2 format
                check_chapter = set_lines[0].strip()
                try:
                    # see if there is a chapter number preceding the set
                    check_chapter = int(check_chapter)
                except:
                    # there was no chapter number, proceed reading a set
                    pass
                else:
                    # there was a chapter number so check that the chapter number
                    # is the same as the first set in this rate file
                    try:
                        assert check_chapter == self.chapter
                    except:
                        print(f'ERROR: read chapter {check_chapter}, expected chapter {self.chapter} for this rate set.')
                        raise
                    else:
                        # get rid of chapter number so we can read a rate set
                        set_lines.pop(0)

                # sets are 3 lines long
                s1 = set_lines.pop(0)
                s2 = set_lines.pop(0)
                s3 = set_lines.pop(0)

                # first line of a set has up to 6 nuclei, then the label,
                # and finally the Q value

                # get rid of first 5 spaces
                s1 = s1[5:]

                # next follows 6 fields of 5 characters containing nuclei
                # the 6 fields are padded with spaces
                f = []
                for i in range(6):
                    ni = s1[:5]
                    s1 = s1[5:]
                    ni = ni.strip()
                    if ni:
                        f.append(ni)

                # next come 8 spaces, so get rid of them
                s1 = s1[8:]

                # next is a 4-character set label and 2 character flags
                labelprops = s1[:6]
                s1 = s1[6:]

                # next come 3 spaces
                s1 = s1[3:]

                # next comes a 12 character Q value followed by 10 spaces
                Q = float(s1.strip())

                if first:
                    self.Q = Q

                    try:
                        # what's left are the nuclei -- their interpretation
                        # depends on the chapter
                        if self.chapter == 1:
                            # e1 -> e2
                            self.reactants.append(Nucleus(f[0]))
                            self.products.append(Nucleus(f[1]))

                        elif self.chapter == 2:
                            # e1 -> e2 + e3
                            self.reactants.append(Nucleus(f[0]))
                            self.products += [Nucleus(f[1]), Nucleus(f[2])]

                        elif self.chapter == 3:
                            # e1 -> e2 + e3 + e4
                            self.reactants.append(Nucleus(f[0]))
                            self.products += [Nucleus(f[1]), Nucleus(f[2]), Nucleus(f[3])]

                        elif self.chapter == 4:
                            # e1 + e2 -> e3
                            self.reactants += [Nucleus(f[0]), Nucleus(f[1])]
                            self.products.append(Nucleus(f[2]))

                        elif self.chapter == 5:
                            # e1 + e2 -> e3 + e4
                            self.reactants += [Nucleus(f[0]), Nucleus(f[1])]
                            self.products += [Nucleus(f[2]), Nucleus(f[3])]

                        elif self.chapter == 6:
                            # e1 + e2 -> e3 + e4 + e5
                            self.reactants += [Nucleus(f[0]), Nucleus(f[1])]
                            self.products += [Nucleus(f[2]), Nucleus(f[3]), Nucleus(f[4])]

                        elif self.chapter == 7:
                            # e1 + e2 -> e3 + e4 + e5 + e6
                            self.reactants += [Nucleus(f[0]), Nucleus(f[1])]
                            self.products += [Nucleus(f[2]), Nucleus(f[3]),
                                              Nucleus(f[4]), Nucleus(f[5])]

                        elif self.chapter == 8:
                            # e1 + e2 + e3 -> e4
                            self.reactants += [Nucleus(f[0]), Nucleus(f[1]), Nucleus(f[2])]
                            self.products.append(Nucleus(f[3]))

                        elif self.chapter == 9:
                            # e1 + e2 + e3 -> e4 + e5
                            self.reactants += [Nucleus(f[0]), Nucleus(f[1]), Nucleus(f[2])]
                            self.products += [Nucleus(f[3]), Nucleus(f[4])]

                        elif self.chapter == 10:
                            # e1 + e2 + e3 + e4 -> e5 + e6
                            self.reactants += [Nucleus(f[0]), Nucleus(f[1]),
                                               Nucleus(f[2]), Nucleus(f[3])]
                            self.products += [Nucleus(f[4]), Nucleus(f[5])]

                        elif self.chapter == 11:
                            # e1 -> e2 + e3 + e4 + e5
                            self.reactants.append(Nucleus(f[0]))
                            self.products += [Nucleus(f[1]), Nucleus(f[2]),
                                              Nucleus(f[3]), Nucleus(f[4])]
                        else:
                            print(f'Chapter could not be identified in {self.original_source}')
                            assert isinstance(self.chapter, int) and self.chapter <= 11
                    except:
                        # print('Error parsing Rate from {}'.format(self.original_source))
                        raise

                    first = 0

                # the second line contains the first 4 coefficients
                # the third lines contains the final 3
                # we can't just use split() here, since the fields run into one another
                n = 13  # length of the field
                a = [s2[i:i+n] for i in range(0, len(s2), n)]
                a += [s3[i:i+n] for i in range(0, len(s3), n)]

                a = [float(e) for e in a if not e.strip() == ""]
                self.sets.append(SingleSet(a, labelprops=labelprops))
                self._set_label_properties(labelprops)

    def _set_rhs_properties(self):
        """ compute statistical prefactor and density exponent from the reactants. """
        self.prefactor = 1.0  # this is 1/2 for rates like a + a (double counting)
        self.inv_prefactor = 1
        for r in set(self.reactants):
            self.inv_prefactor = self.inv_prefactor * np.math.factorial(self.reactants.count(r))
        self.prefactor = self.prefactor/float(self.inv_prefactor)
        self.dens_exp = len(self.reactants)-1
        if (self.weak_type == 'electron_capture' and not self.tabular):
            self.dens_exp = self.dens_exp + 1

    def _set_screening(self):
        """ determine if this rate is eligible for screening and the nuclei to use. """
        # Tells if this rate is eligible for screening
        # using screenz.f90 provided by StarKiller Microphysics.
        # If not eligible for screening, set to None
        # If eligible for screening, then
        # Rate.ion_screen is a 2-element (3 for 3-alpha) list of Nucleus objects for screening
        self.ion_screen = []
        nucz = [q for q in self.reactants if q.Z != 0]
        if len(nucz) > 1:
            nucz.sort(key=lambda x: x.Z)
            self.ion_screen = []
            self.ion_screen.append(nucz[0])
            self.ion_screen.append(nucz[1])
            if len(nucz) == 3:
                self.ion_screen.append(nucz[2])

        # if the rate is a reverse rate, via detailed balance, then we
        # might actually want to compute the screening based on the
        # reactants of the forward rate that was used in the detailed
        # balance.  Rate.symmetric_screen is what should be used in
        # the screening in this case
        self.symmetric_screen = []
        if self.reverse:
            nucz = [q for q in self.products if q.Z != 0]
            if len(nucz) > 1:
                nucz.sort(key=lambda x: x.Z)
                self.symmetric_screen = []
                self.symmetric_screen.append(nucz[0])
                self.symmetric_screen.append(nucz[1])
                if len(nucz) == 3:
                    self.symmetric_screen.append(nucz[2])
        else:
            self.symmetric_screen = self.ion_screen

    def _set_print_representation(self):
        """ compose the string representations of this Rate. """
        self.string = ""
        self.pretty_string = r"$"

        # put p, n, and alpha second
        treactants = []
        for n in self.reactants:
            if n.raw not in ["p", "he4", "n"]:
                treactants.insert(0, n)
            else:
                treactants.append(n)

        for n, r in enumerate(treactants):
            self.string += f"{r}"
            self.pretty_string += fr"{r.pretty}"
            if not n == len(self.reactants)-1:
                self.string += " + "
                self.pretty_string += r" + "

        self.string += " --> "
        self.pretty_string += r" \rightarrow "

        for n, p in enumerate(self.products):
            self.string += f"{p}"
            self.pretty_string += fr"{p.pretty}"
            if not n == len(self.products)-1:
                self.string += " + "
                self.pretty_string += r" + "

        self.pretty_string += r"$"

        if not self.fname:
            # This is used to determine which rates to detect as the same reaction
            # from multiple sources in a Library file, so it should not be unique
            # to a given source, e.g. wc12, but only unique to the reaction.
            reactants_str = '_'.join([repr(nuc) for nuc in self.reactants])
            products_str = '_'.join([repr(nuc) for nuc in self.products])
            self.fname = f'{reactants_str}__{products_str}'
            if self.weak:
                self.fname = self.fname + f'__weak__{self.weak_type}'

    def get_rate_id(self):
        """ Get an identifying string for this rate.
        Don't include resonance state since we combine resonant and
        non-resonant versions of reactions. """

        srev = ''
        if self.reverse:
            srev = 'reverse'

        sweak = ''
        if self.weak:
            sweak = 'weak'

        ssrc = 'reaclib'
        if self.tabular:
            ssrc = 'tabular'

        return f'{self.__repr__()} <{self.label.strip()}_{ssrc}_{sweak}_{srev}>'

    def heaviest(self):
        """
        Return the heaviest nuclide in this Rate.

        If two nuclei are tied in mass number, return the one with the
        lowest atomic number.
        """
        nuc = self.reactants[0]
        for n in self.reactants + self.products:
            if n.A > nuc.A or (n.A == nuc.A and n.Z < nuc.Z):
                nuc = n
        return nuc

    def lightest(self):
        """
        Return the lightest nuclide in this Rate.

        If two nuclei are tied in mass number, return the one with the
        highest atomic number.
        """
        nuc = self.reactants[0]
        for n in self.reactants + self.products:
            if n.A < nuc.A or (n.A == nuc.A and n.Z > nuc.Z):
                nuc = n
        return nuc

    def get_tabular_rate(self):
        """read the rate data from .dat file """

        # find .dat file and read it
        self.table_path = _find_rate_file(self.table_file)
        tabular_file = open(self.table_path)
        t_data = tabular_file.readlines()
        tabular_file.close()

        # delete header lines
        del t_data[0:self.table_header_lines]

        # change the list ["1.23 3.45 5.67\n"] into the list ["1.23","3.45","5.67"]
        t_data2d = []
        for tt in t_data:
            t_data2d.append(re.split(r"[ ]", tt.strip('\n')))

        # delete all the "" in each element of data1
        for tt2d in t_data2d:
            while '' in tt2d:
                tt2d.remove('')

        while [] in t_data2d:
            t_data2d.remove([])

        self.tabular_data_table = np.array(t_data2d)

    def eval(self, T, rhoY = None):
        """ evauate the reaction rate for temperature T """

        if self.tabular:
            data = self.tabular_data_table.astype(np.float)
            # find the nearest value of T and rhoY in the data table
            T_nearest = (data[:,1])[np.abs((data[:,1]) - T).argmin()]
            rhoY_nearest = (data[:,0])[np.abs((data[:,0]) - rhoY).argmin()]
            inde = np.where((data[:,1]==T_nearest)&(data[:,0]==rhoY_nearest))[0][0]
            r = data[inde][5]

        else:
            tf = Tfactors(T)
            r = 0.0
            for s in self.sets:
                f = s.f()
                r += f(tf)

        return r

    def get_rate_exponent(self, T0):
        """
        for a rate written as a power law, r = r_0 (T/T0)**nu, return
        nu corresponding to T0
        """

        # nu = dln r /dln T, so we need dr/dT
        r1 = self.eval(T0)
        dT = 1.e-8*T0
        r2 = self.eval(T0 + dT)

        drdT = (r2 - r1)/dT
        return (T0/r1)*drdT

    def plot(self, Tmin=1.e8, Tmax=1.6e9, rhoYmin=3.9e8, rhoYmax=2.e9):
        """plot the rate's temperature sensitivity vs temperature"""

        if self.tabular:
            data = self.tabular_data_table.astype(np.float) # convert from str to float

            inde1 = data[:,1]<=Tmax
            inde2 = data[:,1]>=Tmin
            inde3 = data[:,0]<=rhoYmax
            inde4 = data[:,0]>=rhoYmin
            data_heatmap = data[inde1&inde2&inde3&inde4].copy()

            rows, row_pos = np.unique(data_heatmap[:, 0], return_inverse=True)
            cols, col_pos = np.unique(data_heatmap[:, 1], return_inverse=True)
            pivot_table = np.zeros((len(rows), len(cols)), dtype=data_heatmap.dtype)
            try:
                pivot_table[row_pos, col_pos] = np.log10(data_heatmap[:, 5])
            except ValueError:
                print("Divide by zero encountered in log10\nChange the scale of T or rhoY")

            _, ax = plt.subplots(figsize=(10,10))
            im = ax.imshow(pivot_table, cmap='jet')
            plt.colorbar(im)

            plt.xlabel("$T$ [K]")
            plt.ylabel("$\\rho Y$ [g/cm$^3$]")
            ax.set_title(fr"{self.pretty_string}"+
                         "\n"+"electron-capture/beta-decay rate in log10(1/s)")
            ax.set_yticks(range(len(rows)))
            ax.set_yticklabels(rows)
            ax.set_xticks(range(len(cols)))
            ax.set_xticklabels(cols)
            plt.setp(ax.get_xticklabels(), rotation=90, ha="right",rotation_mode="anchor")
            plt.gca().invert_yaxis()
            plt.show()

        else:
            temps = np.logspace(np.log10(Tmin), np.log10(Tmax), 100)
            r = np.zeros_like(temps)

            for n, T in enumerate(temps):
                r[n] = self.eval(T)

            plt.loglog(temps, r)
            plt.xlabel(r"$T$")

            if self.dens_exp == 0:
                plt.ylabel(r"\tau")
            elif self.dens_exp == 1:
                plt.ylabel(r"$N_A <\sigma v>$")
            elif self.dens_exp == 2:
                plt.ylabel(r"$N_A^2 <n_a n_b n_c v>$")

            plt.title(fr"{self.pretty_string}")
            plt.show()

class RatePair:
    """the forward and reverse rates for a single reaction sequence.
    Forward rates are those with Q >= 0.

    :var forward: the forward reaction Rate object
    :var reverse: the reverse reaction Rate object

    """

    def __init__(self, forward=None, reverse=None):
        self.forward = forward
        self.reverse = reverse

    def __repr__(self):
        return f"forward: {self.forward} ; reverse: {self.reverse}"

    def __lt__(self, other):
        if self.forward is not None and other.forward is not None:
            return self.forward < other.forward
        elif self.forward is None:
            return False
        return True

    def __eq__(self, other):
        return self.forward == other.forward and self.reverse == other.reverse
