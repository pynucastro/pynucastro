"""
Classes and methods to interface with files storing rate data.
"""

import os
import re
import io
import numpy as np
import matplotlib.pyplot as plt
import collections

from pynucastro.nucdata import UnidentifiedElement, PeriodicTable


def list_known_rates():
    """ list the rates found in the library """

    lib_path = "{}/../library/".format(os.path.dirname(__file__))

    for _, _, filenames in os.walk(lib_path):
        for f in filenames:
            # skip over files that are not rate files
            if f.endswith(".md") or f.endswith(".dat"):
                continue
            try:
                lib = Library(f)
            except:
                continue
            else:
                print("{:32} : ".format(f))
                for r in lib.get_rates():
                    print("                                 : {}".format(r))

import numba

Tfactor_spec = [
('T9', numba.float64),
('T9i', numba.float64),
('T913', numba.float64),
('T913i', numba.float64),
('T953', numba.float64),
('lnT9', numba.float64)
]

@numba.jitclass(Tfactor_spec)
class Tfactors(object):
    """ precompute temperature factors for speed """

    def __init__(self, T):
        """ return the Tfactors object.  Here, T is temperature in Kelvin """
        self.T9 = T/1.e9
        self.T9i = 1.0/self.T9
        self.T913i = self.T9i**(1./3.)
        self.T913 = self.T9**(1./3.)
        self.T953 = self.T9**(5./3.)
        self.lnT9 = np.log(self.T9)


class SingleSet(object):
    """ a set in Reaclib is one piece of a rate, in the form

        lambda = exp[ a_0 + sum_{i=1}^5  a_i T_9**(2i-5)/3  + a_6 log T_9]

        A single rate in Reaclib can be composed of multiple sets
    """

    def __init__(self, a, labelprops=None):
        """here a is iterable (e.g., list or numpy array), storing the
           coefficients, a0, ..., a6

        """
        self.a = a
        self.labelprops = labelprops
        self._update_label_properties()

    def _update_label_properties(self):
        """ Set label and flags indicating Set is resonant,
            weak, or reverse. """
        assert(type(self.labelprops) == str)
        try:
            assert(len(self.labelprops) == 6)
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
            string = "{} += np.exp( ".format(prefix)
        else:
            string = "{} = np.exp( ".format(prefix)
        string += " {}".format(self.a[0])
        if not self.a[1] == 0.0:
            string += " + {}*tf.T9i".format(self.a[1])
        if not self.a[2] == 0.0:
            string += " + {}*tf.T913i".format(self.a[2])
        if not self.a[3] == 0.0:
            string += " + {}*tf.T913".format(self.a[3])
        if not (self.a[4] == 0.0 and self.a[5] == 0.0 and self.a[6] == 0.0):
            string += "\n{}         ".format(len(prefix)*" ")
        if not self.a[4] == 0.0:
            string += " + {}*tf.T9".format(self.a[4])
        if not self.a[5] == 0.0:
            string += " + {}*tf.T953".format(self.a[5])
        if not self.a[6] == 0.0:
            string += " + {}*tf.lnT9".format(self.a[6])
        string += ")"
        return string


class UnsupportedNucleus(BaseException):
    def __init__(self):
        return


class Nucleus(object):
    """
    a nucleus that participates in a reaction -- we store it in a
    class to hold its properties, define a sorting, and give it a
    pretty printing string

    """
    def __init__(self, name):
        name = name.lower()
        self.raw = name

        # element symbol and atomic weight
        if name == "p":
            self.el = "h"
            self.A = 1
            self.short_spec_name = "h1"
        elif name == "d":
            self.el = "h"
            self.A = 2
            self.short_spec_name = "h2"
        elif name == "t":
            self.el = "h"
            self.A = 3
            self.short_spec_name = "h3"
        elif name == "a":
            #this is a convenience, enabling the use of a commonly-used alias:
            #    He4 --> \alpha --> "a" , e.g. c12(a,g)o16
            self.el ="he"
            self.A = 4
            self.short_spec_name = "he4"
            self.raw = "he4"
        elif name == "n":
            self.el = "n"
            self.A = 1
            self.Z = 0
            self.N = 1
            self.short_spec_name = "n"
            self.spec_name = "neutron"
            self.pretty = r"\mathrm{{{}}}".format(self.el)
        else:
            e = re.match(r"([a-zA-Z]*)(\d*)", name)
            self.el = e.group(1).title()  # chemical symbol
            assert(self.el)
            try:
                self.A = int(e.group(2))
            except:
                if (name.strip() == 'al-6' or
                    name.strip() == 'al*6'):
                    raise UnsupportedNucleus()
                else:
                    raise
            assert(self.A >= 0)
            self.short_spec_name = name

        # use lowercase element abbreviation regardless the case of the input
        self.el = self.el.lower()

        # atomic number comes from periodic table
        if name != "n":
            try:
                i = PeriodicTable.lookup_abbreviation(self.el)
            except UnidentifiedElement:
                print('Could not identify element: {}'.format(self.el))
                raise
            except:
                raise
            else:
                self.Z = i.Z
                assert(type(self.Z) == int)
                assert(self.Z >= 0)
                self.N = self.A - self.Z
                assert(type(self.N) == int)
                assert(self.N >= 0)

                # long name
                self.spec_name = '{}-{}'.format(i.name, self.A)

                # latex formatted style
                self.pretty = r"{{}}^{{{}}}\mathrm{{{}}}".format(self.A, self.el.capitalize())

    def __repr__(self):
        return self.raw

    def __hash__(self):
        return hash(self.__repr__())

    def __eq__(self, other):
        return self.el == other.el and \
               self.Z == other.Z and self.A == other.A

    def __lt__(self, other):
        if not self.Z == other.Z:
            return self.Z < other.Z
        else:
            return self.A < other.A


class Library(object):
    """
    A Library is a Rate container that reads a single file
    containing one or many Reaclib rates, possibly containing multiple
    sets per rate.

    The Library class also implements searching based on rules
    specified by RateFilter objects.
    """

    pynucastro_dir = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
    pynucastro_rates_dir = os.path.join(pynucastro_dir,
                                        'library')
    pynucastro_tabular_dir = os.path.join(pynucastro_rates_dir,
                                          'tabular')

    def __init__(self, libfile=None, rates=None, read_library=True):
        self._library_file = libfile
        if rates:
            self._rates = None
            if isinstance(rates, Rate):
                rates = [rates]
            assert (isinstance(rates, dict) or isinstance(rates, list)), "ERROR: rates in Library constructor must be a Rate object, list of Rate objects, or dictionary of Rate objects keyed by Rate.get_rate_id()"
            if isinstance(rates, dict):
                self._rates = rates
            elif isinstance(rates, list):
                self._add_from_rate_list(rates)
        else:
            self._rates = collections.OrderedDict()
        self._library_source_lines = []

        if self._library_file and read_library:
            self._library_file = self._find_rate_file(self._library_file)
            self._read_library_file()

    def heaviest(self):
        """ Return the heaviest nuclide in this library. """
        nuc = None
        for id, r in self._rates.items():
            rnuc = r.heaviest()
            if nuc:
                if rnuc.A > nuc.A or (rnuc.A == nuc.A and rnuc.Z < nuc.Z):
                    nuc = rnuc
            else:
                nuc = rnuc
        return nuc

    def lightest(self):
        """ Return the lightest nuclide in this library. """
        nuc = None
        for id, r in self._rates.items():
            rnuc = r.lightest()
            if nuc:
                if rnuc.A < nuc.A or (rnuc.A == nuc.A and rnuc.Z > nuc.Z):
                    nuc = rnuc
            else:
                nuc = rnuc
        return nuc

    def _add_from_rate_list(self, ratelist):
        """ Add to the rate dictionary from the supplied list of Rate objects. """
        if not self._rates:
            self._rates = collections.OrderedDict()
        for r in ratelist:
            id = r.get_rate_id()
            assert (not id in self._rates), "ERROR: supplied a Rate object already in the Library."
            self._rates[id] = r

    @classmethod
    def _find_rate_file(self, ratename):
        """locate the Reaclib or tabular rate or library file given its name.  Return
        None if the file cannot be located, otherwise return its path."""

        # check to see if the rate file is in the working dir or
        # is already the full path
        x = ratename
        if os.path.isfile(x):
            return os.path.realpath(x)

        # check to see if the rate file is in pynucastro/library
        x = os.path.join(self.pynucastro_rates_dir, ratename)
        if os.path.isfile(x):
            return os.path.realpath(x)

        # check to see if the rate file is in pynucastro/library/tabular
        x = os.path.join(self.pynucastro_tabular_dir, ratename)
        if os.path.isfile(x):
            return os.path.realpath(x)

        # notify user we can't find the file
        raise Exception('File {} not found in the working directory, {}, or {}'.format(
            ratename, self.pynucastro_rates_dir, self.pynucastro_tabular_dir))

    def _read_library_file(self):
        # loop through library file, read lines
        try:
            flib = open(self._library_file, 'r')
        except:
            print('Could not open file {}'.format(self._library_file))
            raise
        for line in flib:
            ls = line.rstrip('\n')
            if ls.strip():
                self._library_source_lines.append(ls)
        flib.close()

        # identify distinct rates from library lines
        current_chapter = None
        while True:
            if len(self._library_source_lines) == 0:
                break

            # Check to see if there is a chapter ID, if not then use current_chapter
            # (for Reaclib v1 formatted library files)
            line = self._library_source_lines[0].strip()
            chapter = None
            if line == 't' or line == 'T':
                chapter = 't'
                self._library_source_lines.pop(0)
            else:
                try:
                    chapter = int(line)
                except:
                    # we can't interpret line as a chapter so use current_chapter
                    try:
                        assert(current_chapter)
                    except:
                        print('ERROR: malformed library file {}, cannot identify chapter.'.format(self._library_file))
                        raise
                    else:
                        chapter = current_chapter
                else:
                    self._library_source_lines.pop(0)
            current_chapter = chapter

            rlines = None
            if chapter == 't':
                rlines = [self._library_source_lines.pop(0) for i in range(5)]
            elif type(chapter) == int:
                rlines = [self._library_source_lines.pop(0) for i in range(3)]
            if rlines:
                sio = io.StringIO('\n'.join(['{}'.format(chapter)] +
                                            rlines))
                #print(sio.getvalue())
                try:
                    r = Rate(sio, rfile_path=self._library_file)
                except UnsupportedNucleus:
                    pass
                except:
                    raise
                else:
                    id = r.get_rate_id()
                    if id in self._rates:
                        self._rates[id] = self._rates[id] + r
                    else:
                        self._rates[id] = r

    def __repr__(self):
        """ Return a string containing the rates IDs in this library. """
        rstrings = []
        for id, r in self._rates.items():
            rstrings.append('{}    ({})'.format(r, id))
        return '\n'.join(rstrings)

    def __add__(self, other):
        """ Add two libraries to get a library containing rates from both. """
        new_rates = self._rates
        for id, r in other._rates.items():
            try:
                assert id not in new_rates
            except:
                if r != new_rates[id]:
                    print('ERROR: rate {} defined differently in libraries {} and {}\n'.format(r, self._library_file, other._library_file))
                    raise
            else:
                new_rates[id] = r
        new_library = Library(libfile='{} + {}'.format(self._library_file, other._library_file),
                              rates=new_rates,
                              read_library=False)
        return new_library

    def get_rates(self):
        """ Return a list of the rates in this library. """
        rlist = [r for id, r in self._rates.items()]
        return rlist

    def get_rate(self, id):
        """ Return a rate matching the id provided. """
        try:
            return self._rates[id]
        except:
            print("ERROR: rate identifier does not match a rate in this library.")
            raise

    def linking_nuclei(self, nuclist, with_reverse=True):
        """
        Return a Library object containing the rates linking the
        nuclei provided in the list of Nucleus objects or nucleus abbreviations 'nuclist'.

        If with_reverse is True, then include reverse rates. Otherwise
        include only forward rates.
        """

        if type(nuclist) == Nucleus or type(nuclist) == str:
            nuclist = [nuclist]
        else:
            try:
                nuclist = list(nuclist)
            except:
                raise

        nucleus_list = []
        for nuc in nuclist:
            if type(nuc) == Nucleus:
                nucleus_list.append(nuc)
            else:
                try:
                    anuc = Nucleus(nuc)
                except:
                    raise
                else:
                    nucleus_list.append(anuc)

        # Get the set of rates for which any Nucleus in nucleus_list
        # appears as either reactant or product.
        rate_filters = []
        for nuc in nucleus_list:
            rate_filters.append(RateFilter(reactants=nuc, exact=False))
            rate_filters.append(RateFilter(products=nuc, exact=False))
        triage_library = self.filter(rate_filters)

        # Discard any of this set of rates for which nuclei appear not
        # in nucleus_list
        filtered_rates = []
        for r in triage_library.get_rates():
            include = True
            for nuc in r.reactants:
                if nuc not in nucleus_list:
                    include = False
                    break
            for nuc in r.products:
                if nuc not in nucleus_list:
                    include = False
                    break
            if not with_reverse and r.reverse:
                include = False
            if include:
                filtered_rates.append(r)

        # Return library containing the filtered rates
        return Library(rates=filtered_rates)

    def filter(self, filter_spec):
        """
        filter_specs should be an iterable of RateFilter objects or a
        single RateFilter object. Library.filter yields all rates
        matching any RateFilter in filter_specs.  If RateFilter.exact,
        then return rates with exactly the reactants or products
        passed in as arguments.
        """
        if type(filter_spec) == RateFilter:
            filter_specifications = [filter_spec]
        else:
            try:
                iter(filter_spec)
            except:
                raise
            else:
                filter_specifications = filter_spec
        matching_rates = collections.OrderedDict()
        for id, r in self._rates.items():
            for f in filter_specifications:
                if f.matches(r):
                    matching_rates[id] = r
                    break
        if matching_rates:
            return Library(libfile=self._library_file,
                           rates=matching_rates,
                           read_library=False)
        else:
            return None


class RateFilter(object):
    """RateFilter filters out a specified rate or set of rates
    
    A RateFilter stores selection rules specifying a rate or group of
    rates to assist in searching for rates stored in a Library.
    """

    def __init__(self, reactants=None, products=None, exact=True,
                 reverse=None, min_reactants=None, max_reactants=None,
                 min_products=None, max_products=None):
        """Create a new RateFilter with the given selection rules

        Keyword Arguments:
            reactants -- Description of the reactants as one of:
                1. a list of Nucleus objects
                2. a list of string descriptions of reactant nuclides
                   these strings must be parsable by Nucleus
                3. a single reactant Nucleus
                4. a single string description of the reactant nuclide
            products  -- Description of the products in same form as above
            exact     -- boolean, 
                         if True, products or reactants must match exactly [default]
                         if False, then all products or reactants must be found
                         in a comparison rate, but the comparison may contain
                         additional products or reactants
            reverse   -- boolean,
                         if True, only match reverse-derived rates
                         if False, only match directly-derived rates
                         if None, you don't care, match both [default]
            min_reactants -- int, match Rates that have at least this many reactants
            min_products  -- int, match Rates that have at least this many products
            max_reactants -- int, match Rates that have no more than this many reactants
            max_products  -- int, match Rates that have no more than this many products
        
        Examples:
            Create a filter that finds all proton capture and proton-burning reactions
            in a Library instance my_library::
                >>> pcap_filter = RateFilter(reactants='p', exact=False)
                >>> pcap_library = my_library.filter(pcap_filter)
            or you can use Nucleus::
                >>> pcap_filter = RateFilter(reactants=Nucleus('p'), exact=False)
                >>> pcap_library = my_library.filter(pcap_filter)

            Create a filter that finds C12 (a,g) O16 
            Notes:
                + photons/gammas are not treated as nuclides, so they cannot be
                a reactant or product
                + this rate is in the ReacLib library used here as 
                O16 --> He4 C12 -- you need to know how your library treats rates::
                    >>> cago_filter = RateFilter(reactants='o16', products=['c12', 'a'])
                    >>> cago_library = my_library.filter(cago_filter)
        """
        self.reactants = []
        self.products = []
        self.exact = exact
        self.reverse = reverse
        self.min_reactants = min_reactants
        self.min_products = min_products
        self.max_reactants = max_reactants
        self.max_products = max_products

        if reactants:
            if type(reactants) == Nucleus or type(reactants) == str:
                reactants = [reactants]
            self.reactants = [self._cast_nucleus(r) for r in reactants]
        if products:
            if type(products) == Nucleus or type(products) == str:
                products = [products]
            self.products = [self._cast_nucleus(r) for r in products]

    @staticmethod
    def _cast_nucleus(r):
        """ Make sure r is of type Nucleus. """
        if not type(r) == Nucleus:
            try:
                rnuc = Nucleus(r)
            except:
                raise
            else:
                return rnuc
        else:
            return r

    @staticmethod
    def _contents_equal(a, b):
        """
        Return True if the contents of a and b exactly match, ignoring ordering.
        If either a or b is None, return True only if both a and b are None.
        """
        if a and b:
            return collections.Counter(a) == collections.Counter(b)
        else:
            return (not a) and (not b)

    @staticmethod
    def _compare_nuclides(test, reference, exact=True):
        """
        test and reference should be iterables of Nucleus objects.
        If an exact match is desired, test and reference should exactly match, ignoring ordering.
        Otherwise, return True only if every element of test appears at least one time in reference.
        """
        matches = True
        if exact:
            matches = RateFilter._contents_equal(test, reference)
        else:
            for nuc in test:
                if not (nuc in reference):
                    matches = False
                    break
        return matches

    def matches(self, r):
        """ Given a Rate r, see if it matches this RateFilter. """
        matches_reactants = True
        matches_products = True
        matches_reverse = True
        matches_min_reactants = True
        matches_min_products = True
        matches_max_reactants = True
        matches_max_products = True
        if self.reactants:
            matches_reactants = self._compare_nuclides(self.reactants, r.reactants, self.exact)
        if self.products:
            matches_products = self._compare_nuclides(self.products, r.products, self.exact)
        if type(self.reverse) == type(True):
            matches_reverse = self.reverse == r.reverse
        if type(self.min_reactants) == int:
            matches_min_reactants = len(r.reactants) >= self.min_reactants
        if type(self.min_products) == int:
            matches_min_products = len(r.products) >= self.min_products
        if type(self.max_reactants) == int:
            matches_max_reactants = len(r.reactants) <= self.max_reactants
        if type(self.max_products) == int:
            matches_max_products = len(r.products) <= self.max_products
        return (matches_reactants and matches_products and matches_reverse and
                matches_min_reactants and matches_max_reactants and
                matches_min_products and matches_max_products)

    def invert(self):
        """ Return a RateFilter matching the inverse rate. """
        newfilter = RateFilter(reactants=self.products,
                               products=self.reactants,
                               exact=self.exact,
                               reverse=self.reverse,
                               min_reactants=self.min_products,
                               max_reactants=self.max_products,
                               min_products=self.min_reactants,
                               max_products=self.max_reactants)
        return newfilter


class Rate(object):
    """ a single Reaclib rate, which can be composed of multiple sets """
    def __init__(self, rfile=None, rfile_path=None, chapter=None, original_source=None,
                 reactants=None, products=None, sets=None, labelprops=None, Q=None):
        """ rfile can be either a string specifying the path to a rate file or
        an io.StringIO object from which to read rate information. """

        self.rfile_path = rfile_path
        self.rfile = None

        if type(rfile) == str:
            self.rfile_path = Library._find_rate_file(rfile)
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

        self.Q = Q

        if type(rfile) == str:
            # read in the file, parse the different sets and store them as
            # SingleSet objects in sets[]
            f = open(self.rfile_path, "r")
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

    def __add__(self, other):
        """Combine the sets of two Rate objects if they describe the same
           reaction. Must be Reaclib rates."""
        assert(self.reactants == other.reactants)
        assert(self.products == other.products)
        assert(self.chapter == other.chapter)
        assert(type(self.chapter) == int)
        assert(self.label == other.label)
        assert(self.weak == other.weak)
        assert(self.weak_type == other.weak_type)
        assert(self.tabular == other.tabular)
        assert(self.reverse == other.reverse)

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
        assert(type(self.labelprops) == str)
        try:
            assert(len(self.labelprops) == 6)
        except:
            assert(self.labelprops == 'tabular')
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
                print('Nucleus objects not be identified in {}'.format(self.original_source))
                raise

            self.table_file = s2.strip()
            self.table_header_lines = int(s3.strip())
            self.table_rhoy_lines = int(s4.strip())
            self.table_temp_lines = int(s5.strip())
            self.table_num_vars = 6 # Hard-coded number of variables in tables for now.
            self.table_index_name = 'j_{}_{}'.format(self.reactants[0], self.products[0])
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
                        assert(check_chapter == self.chapter)
                    except:
                        print('ERROR: read chapter {}, expected chapter {} for this rate set.'.format(check_chapter, self.chapter))
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
                Q = s1.strip()

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
                            print('Chapter could not be identified in {}'.format(self.original_source))
                            assert(type(self.chapter) == int and self.chapter <= 11)
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
        # Rate.ion_screen is a 2-element list of Nucleus objects for screening
        self.ion_screen = []
        nucz = []
        for parent in self.reactants:
            if parent.Z != 0:
                nucz.append(parent)
        if len(nucz) > 1:
            nucz.sort(key=lambda x: x.Z)
            self.ion_screen = []
            self.ion_screen.append(nucz[0])
            self.ion_screen.append(nucz[1])

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
            self.string += "{}".format(r)
            self.pretty_string += r"{}".format(r.pretty)
            if not n == len(self.reactants)-1:
                self.string += " + "
                self.pretty_string += r" + "

        self.string += " --> "
        self.pretty_string += r" \rightarrow "

        for n, p in enumerate(self.products):
            self.string += "{}".format(p)
            self.pretty_string += r"{}".format(p.pretty)
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
            self.fname = '{}__{}'.format(reactants_str, products_str)
            if self.weak:
                self.fname = self.fname + '__weak__{}'.format(self.weak_type)

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

        return '{} <{}_{}_{}_{}>'.format(self.__repr__(), self.label.strip(),
                                         ssrc, sweak, srev)

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
        self.table_path = Library._find_rate_file(self.table_file)
        tabular_file = open(self.table_path,"r")
        t_data = tabular_file.readlines()  
        tabular_file.close()
        
        # delete header lines
        del t_data[0:self.table_header_lines]  
        
        # change the list ["1.23 3.45 5.67\n"] into the list ["1.23","3.45","5.67"]
        t_data2d = []
        for i in range(len(t_data)):
            t_data2d.append(re.split(r"[ ]",t_data[i].strip('\n')))
        
        # delete all the "" in each element of data1
        for i in range(len(t_data2d)):
            while '' in t_data2d[i]:
                t_data2d[i].remove('')

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
                plot("Divide by zero encountered in log10\nChange the scale of T or rhoY")
            
            fig, ax = plt.subplots(figsize=(10,10))
            im = ax.imshow(pivot_table, cmap='jet')
            plt.colorbar(im)
            
            plt.xlabel("$T$ [K]")
            plt.ylabel("$\\rho Y$ [g/cm$^3$]")
            ax.set_title(r"{}".format(self.pretty_string)+
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
                
            plt.title(r"{}".format(self.pretty_string))
            plt.show()
            
        
