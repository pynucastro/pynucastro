import collections
import io
import os
import re

from pynucastro.nucdata import Nucleus, UnsupportedNucleus
from pynucastro.rates.rate import (DerivedRate, Rate, ReacLibRate, TabularRate,
                                   _find_rate_file, load_rate)


def list_known_rates():
    """ list the rates found in the library """

    lib_path = f"{os.path.dirname(__file__)}/../library/"

    for _, _, filenames in os.walk(lib_path):
        for f in filenames:
            # skip over files that are not rate files
            if f.endswith(".md") or f.endswith(".dat") or f.endswith(".py") or f.endswith(".ipynb"):
                continue
            try:
                lib = Library(f)
            except Exception:  # pylint: disable=broad-except
                continue
            print(f"{f:32} : ")
            for r in lib.get_rates():
                print(f"                                 : {r}")


def _rate_name_to_nuc(name):

    # first try to interpret name as A(x,y)B
    rate_str = re.compile(r"([A-Za-z0-9]+)\(([A-Za-z0-9_]*),([A-Za-z0-9_]*)\)([A-Za-z0-9]+)",
                          re.IGNORECASE)
    nucs = rate_str.search(name)
    try:
        _r = [nucs.group(1), nucs.group(2)]
        _p = [nucs.group(3), nucs.group(4)]
    except AttributeError:
        return None

    # now try to make nuclei objects.
    reactants = []
    for nuc in _r:
        if nuc == "":
            continue
        try:
            n = Nucleus(nuc)
            reactants.append(n)
        except (ValueError, AssertionError):
            # we need to interpret some things specially
            if nuc.lower() in ["e", "nu", "_", "g", "gamma"]:
                # first electrons and neutrins, and nothing
                continue
            elif nuc.lower() == "aa":
                reactants.append(Nucleus("he4"))
                reactants.append(Nucleus("he4"))
            else:
                raise

    products = []
    for nuc in _p:
        if nuc == "":
            continue
        try:
            n = Nucleus(nuc)
            products.append(n)
        except (ValueError, AssertionError):
            # we need to interpret some things specially
            if nuc.lower() in ["e", "nu", "_", "g", "gamma"]:
                # first electrons and neutrinos, gammas, and nothing
                continue
            elif nuc.lower() == "aa":
                products.append(Nucleus("he4"))
                products.append(Nucleus("he4"))
            else:
                raise

    return reactants, products


class Library:
    """
    A Library is a Rate container that reads a single file
    containing one or many Reaclib rates, possibly containing multiple
    sets per rate.

    The Library class also implements searching based on rules
    specified by RateFilter objects.
    """

    def __init__(self, libfile=None, rates=None):
        self._library_file = libfile
        if rates:
            self._rates = None
            if isinstance(rates, Rate):
                rates = [rates]
            if isinstance(rates, dict):
                self._rates = rates
            elif isinstance(rates, (list, set)):
                self._add_from_rate_list(rates)
            else:
                raise TypeError("rates in Library constructor must be a Rate object, list of Rate objects, or dictionary of Rate objects keyed by Rate.get_rate_id()")
        else:
            self._rates = {}
        self._library_source_lines = collections.deque()

        if self._library_file:
            self._library_file = _find_rate_file(self._library_file)
            self._read_library_file()

    def heaviest(self):
        """ Return the heaviest nuclide in this library. """
        nuc = None
        for _, r in self._rates.items():
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
        for _, r in self._rates.items():
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
            self._rates = {}
        for r in ratelist:
            rid = r.get_rate_id()
            if rid in self._rates:
                raise ValueError(f"supplied a Rate object already in the Library: {r}")
            self._rates[rid] = r

    def _read_library_file(self):
        # loop through library file, read lines
        with open(self._library_file) as flib:
            for line in flib:
                ls = line.rstrip('\n')
                if ls.strip():
                    self._library_source_lines.append(ls)

        # identify distinct rates from library lines
        current_chapter = None
        while True:
            if len(self._library_source_lines) == 0:
                break

            # Check to see if there is a chapter ID, if not then use current_chapter
            # (for Reaclib v1 formatted library files)
            line = self._library_source_lines[0].strip()
            chapter = None
            if line in ('t', 'T'):
                chapter = 't'
                self._library_source_lines.popleft()
            else:
                try:
                    chapter = int(line)
                except (TypeError, ValueError):
                    # we can't interpret line as a chapter so use current_chapter
                    assert current_chapter, f'malformed library file {self._library_file}, cannot identify chapter.'
                    chapter = current_chapter
                else:
                    self._library_source_lines.popleft()
            current_chapter = chapter

            rlines = None
            rate_type = None
            if chapter == 't':
                rlines = [self._library_source_lines.popleft() for i in range(5)]
                rate_type = "tabular"
            elif isinstance(chapter, int):
                rlines = [self._library_source_lines.popleft() for i in range(3)]
                rate_type = "reaclib"
            if rlines:
                sio = io.StringIO('\n'.join([f'{chapter}'] +
                                            rlines))
                #print(sio.getvalue())
                try:
                    if rate_type == "reaclib":
                        r = ReacLibRate(rfile=sio, rfile_path=self._library_file)
                    elif rate_type == "tabular":
                        r = TabularRate(rfile=sio, rfile_path=self._library_file)
                    else:
                        raise NotImplementedError("rate not implemented")
                except UnsupportedNucleus:
                    pass
                else:
                    rid = r.get_rate_id()
                    if rid in self._rates:
                        self._rates[rid] = self._rates[rid] + r
                    else:
                        self._rates[rid] = r

    def __repr__(self):
        """ Return a string containing the rates IDs in this library. """
        rstrings = []
        tmp_rates = [v for k, v in self._rates.items()]
        for r in sorted(tmp_rates):
            if r.Q is not None and r.Q >= 0:
                rstrings.append(f'{r.__repr__():30} [Q = {float(r.Q):6.2f} MeV] ({r.get_rate_id()})')
        for r in sorted(tmp_rates):
            if r.Q is not None and r.Q < 0:
                rstrings.append(f'{r.__repr__():30} [Q = {float(r.Q):6.2f} MeV] ({r.get_rate_id()})')

        for r in sorted(tmp_rates):
            if r.Q is None:
                rstrings.append(f'{r.__repr__():30} ({r.get_rate_id()})')

        return '\n'.join(rstrings)

    def __add__(self, other):
        """ Add two libraries to get a library containing rates from both. """
        new_rates = self._rates
        for rid, r in other._rates.items():
            if rid in new_rates:
                if r != new_rates[rid]:
                    raise ValueError(f'rate {r} defined differently in libraries {self._library_file} and {other._library_file}')
            else:
                new_rates[rid] = r
        new_library = Library(rates=new_rates)
        return new_library

    def __sub__(self, other):
        return self.diff(other)

    def get_num_rates(self):
        """Return the number of rates known to the library"""
        return len(self._rates)

    def get_rates(self):
        """ Return a list of the rates in this library. """
        rlist = [r for _, r in self._rates.items()]
        return rlist

    def get_rate(self, rid):
        """ Return a rate matching the id provided. """
        try:
            return self._rates[rid]
        except KeyError:
            pass

        # fallback to the rate fname
        try:
            r = [q for q in self.get_rates() if q.fname == rid][0]
            return r
        except IndexError:
            raise LookupError(f"rate identifier {rid!r} does not match a rate in this library.") from None

    def get_rate_by_nuclei(self, reactants, products):
        """given a list of reactants and products, return any matching rates"""
        _tmp = [r for r in self.get_rates() if
                sorted(r.reactants) == sorted(reactants) and
                sorted(r.products) == sorted(products)]

        if not _tmp:
            return None
        if len(_tmp) == 1:
            return _tmp[0]
        return _tmp

    def get_rate_by_name(self, name):
        """Given a string representing a rate in the form 'A(x,y)B'
        (or a list of strings for multiple rates) return the Rate
        objects that match from the Library.  If there are multiple
        inputs, then a list of Rate objects is returned.

        """

        if isinstance(name, str):
            rate_name_list = [name]
        else:
            rate_name_list = name

        rates_out = []

        for rname in rate_name_list:
            reactants, products = _rate_name_to_nuc(rname)

            rf = RateFilter(reactants=reactants, products=products)
            _lib = self.filter(rf)
            if _lib is None:
                return None
            rates_out += _lib.get_rates()

        if (len(rates_out)) == 1:
            return rates_out[0]
        return rates_out

    def get_nuclei(self):
        """get the list of unique nuclei"""
        return {nuc for r in self.get_rates() for nuc in r.reactants + r.products}

    def diff(self, other_library):
        """Return a Library containing the rates in this library that are not
        contained in other_library"""

        diff_rates = set(self.get_rates()) - set(other_library.get_rates())
        new_library = Library(rates=diff_rates)
        return new_library

    def remove_rate(self, rate):
        """Manually remove a rate from the library by supplying the
        short name "A(x,y)B, a Rate object, or the rate id"""

        if isinstance(rate, Rate):
            rid = rate.get_rate_id()
            self._rates.pop(rid)
        elif isinstance(rate, str):
            rid = self.get_rate_by_name(rate).get_rate_id()
            self._rates.pop(rid)
        else:
            # we assume that a rate id as provided
            self._rates.pop(rate)

    def add_rate(self, rate):
        """Manually add a rate by giving a Rate object"""

        if isinstance(rate, Rate):
            self._rates[rate.get_rate_id()] = rate
        else:
            raise TypeError("invalid Rate object")

    def linking_nuclei(self, nuclist, with_reverse=True):
        """
        Return a Library object containing the rates linking the
        nuclei provided in the list of Nucleus objects or nucleus abbreviations 'nuclist'.

        If with_reverse is True, then include reverse rates. Otherwise
        include only forward rates.
        """

        if isinstance(nuclist, (Nucleus, str)):
            nuclist = [nuclist]

        nucleus_set = set()
        for nuc in nuclist:
            if isinstance(nuc, Nucleus):
                nucleus_set.add(nuc)
            else:
                anuc = Nucleus(nuc)
                nucleus_set.add(anuc)

        # Discard rates with nuclei that are not in nucleus_set
        filtered_rates = []
        for r in self.get_rates():
            include = True
            for nuc in r.reactants:
                if nuc not in nucleus_set:
                    include = False
                    break
            for nuc in r.products:
                if nuc not in nucleus_set:
                    include = False
                    break
            if not with_reverse and r.reverse:
                include = False
            if include:
                filtered_rates.append(r)

        # create a new library containing the filtered rates
        new_lib = Library(rates=filtered_rates)

        # print out a warning if one of the input nuclei is not linked
        lib_nuclei = new_lib.get_nuclei()
        for nuc in nucleus_set:
            if nuc not in lib_nuclei:
                print(f"warning {nuc} was not able to be linked")

        return new_lib

    def filter(self, filter_spec):
        """
        filter_specs should be an iterable of RateFilter objects or a
        single RateFilter object. Library.filter yields all rates
        matching any RateFilter in filter_specs.  If RateFilter.exact,
        then return rates with exactly the reactants or products
        passed in as arguments.
        """
        if isinstance(filter_spec, RateFilter):
            filter_specifications = [filter_spec]
        else:
            filter_specifications = list(filter_spec)
        matching_rates = {}
        for rid, r in self._rates.items():
            for f in filter_specifications:
                if f.matches(r):
                    matching_rates[rid] = r
                    break
        if matching_rates:
            return Library(rates=matching_rates)
        return None

    def forward(self):
        """
        Select only the forward rates, discarding the inverse rates obtained
        by detailed balance.
        """

        only_fwd_filter = RateFilter(reverse=False)
        only_fwd = self.filter(only_fwd_filter)
        return only_fwd

    def backward(self):
        """
        Select only the reverse rates, obtained by detailed balance.
        """

        only_bwd_filter = RateFilter(reverse=True)
        only_bwd = self.filter(only_bwd_filter)
        return only_bwd

    def derived_forward(self):
        """
        In this library, We exclude the weak and tabular rates from the .foward() library which includes all
        the ReacLib forward reactions.

        In a future PR, we will classify forward reactions as exothermic (Q>0), and reverse by endothermic (Q<0).
        However, ReacLib does not follow this path. If a reaction is measured experimentally (independent of Q),
        they use detailed balance to get the opposite direction. Eventually, I want to classify forward and reverse
        by positive Q and negative Q; however, for testing purposes, making this classification may eventually lead to
        computing the detailed balance twice.

        The idea of derived_forward is to eliminate the reverse and weak, and see if our job gives the same Reaclib
        predictions, checking the NSE convergence with the pf functions. In the future, I want to move this function
        in a unit test.
        """

        collect_rates = []
        onlyfwd = self.forward()

        for r in onlyfwd.get_rates():

            try:
                DerivedRate(rate=r, compute_Q=False, use_pf=False)
            except ValueError:
                continue
            else:
                collect_rates.append(r)

        list1 = Library(rates=collect_rates)
        return list1

    def derived_backward(self, compute_Q=False, use_pf=False):
        """
        This library contains the detailed balance reverse reactions over the selected .derived_forward(),
        computed by hand.
        """

        derived_rates = []
        onlyfwd = self.derived_forward()

        for r in onlyfwd.get_rates():
            try:
                i = DerivedRate(rate=r, compute_Q=compute_Q, use_pf=use_pf)
            except ValueError:
                continue
            else:
                derived_rates.append(i)

        onlybwd = Library(rates=derived_rates)
        return onlybwd


class RateFilter:
    """RateFilter filters out a specified rate or set of rates
    A RateFilter stores selection rules specifying a rate or group of
    rates to assist in searching for rates stored in a Library.
    """

    def __init__(self, reactants=None, products=None, exact=True,
                 reverse=None, min_reactants=None, max_reactants=None,
                 min_products=None, max_products=None, filter_function=None):
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
            filter_function -- callable (Rate -> bool),
                               a callable that can take a single rate as an argument
                               may be used to specify additional criteria, returning
                               True if the rate meets all of them, False otherwise

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
        self.filter_function = filter_function

        if reactants:
            if isinstance(reactants, (Nucleus, str)):
                reactants = [reactants]
            self.reactants = [self._cast_nucleus(r) for r in reactants]
        if products:
            if isinstance(products, (Nucleus, str)):
                products = [products]
            self.products = [self._cast_nucleus(r) for r in products]

    @staticmethod
    def _cast_nucleus(r):
        """ Make sure r is of type Nucleus. """
        if not isinstance(r, Nucleus):
            rnuc = Nucleus(r)
            return rnuc
        return r

    @staticmethod
    def _contents_equal(a, b):
        """
        Return True if the contents of a and b exactly match, ignoring ordering.
        If either a or b is None, return True only if both a and b are None.
        """
        if a and b:
            return len(a) == len(b) and sorted(a) == sorted(b)
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
                if nuc not in reference:
                    matches = False
                    break
        return matches

    def matches(self, r):
        """ Given a Rate r, see if it matches this RateFilter. """
        # do cheaper checks first
        matches_reverse = True
        matches_min_reactants = True
        matches_min_products = True
        matches_max_reactants = True
        matches_max_products = True
        if isinstance(self.reverse, bool):
            matches_reverse = self.reverse == r.reverse
        if isinstance(self.min_reactants, int):
            matches_min_reactants = len(r.reactants) >= self.min_reactants
        if isinstance(self.min_products, int):
            matches_min_products = len(r.products) >= self.min_products
        if isinstance(self.max_reactants, int):
            matches_max_reactants = len(r.reactants) <= self.max_reactants
        if isinstance(self.max_products, int):
            matches_max_products = len(r.products) <= self.max_products
        # exit early if any of these checks failed
        if not (matches_reverse and matches_min_reactants and
                matches_min_products and matches_max_reactants and
                matches_max_products):
            return False
        # now do more expensive checks, and exit immediately if any fail
        if self.reactants:
            if not self._compare_nuclides(self.reactants, r.reactants, self.exact):
                return False
        if self.products:
            if not self._compare_nuclides(self.products, r.products, self.exact):
                return False
        if self.filter_function is not None:
            if not self.filter_function(r):
                return False
        return True

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


class ReacLibLibrary(Library):
    """Load the latest stored version of the ReacLib library and
    return a Library"""

    def __init__(self):
        libfile = 'reaclib_default2_20220329'
        Library.__init__(self, libfile=libfile)


class TabularLibrary(Library):
    """Load all of the tabular rates known and return a Library"""

    def __init__(self):
        # find all of the tabular rates that pynucastro knows about
        # we'll assume that these are of the form *-toki

        lib_path = f"{os.path.dirname(__file__)}/../library/"

        trates = []

        for _, _, filenames in os.walk(lib_path):
            for f in filenames:
                if f.endswith("-toki"):
                    trates.append(load_rate(f))

        Library.__init__(self, rates=trates)
