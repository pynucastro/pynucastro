import os
import io
import collections

from rate import Rate, UnsupportedNucleus, _find_rate_file


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
            except:
                continue
            else:
                print(f"{f:32} : ")
                for r in lib.get_rates():
                    print(f"                                 : {r}")

class Library:
    """
    A Library is a Rate container that reads a single file
    containing one or many Reaclib rates, possibly containing multiple
    sets per rate.

    The Library class also implements searching based on rules
    specified by RateFilter objects.
    """

    def __init__(self, libfile=None, rates=None, read_library=True):
        self._library_file = libfile
        if rates:
            self._rates = None
            if isinstance(rates, Rate):
                rates = [rates]
            assert isinstance(rates, dict) or isinstance(rates, list) or isinstance(rates, set), "ERROR: rates in Library constructor must be a Rate object, list of Rate objects, or dictionary of Rate objects keyed by Rate.get_rate_id()"
            if isinstance(rates, dict):
                self._rates = rates
            elif isinstance(rates, list) or isinstance(rates, set):
                self._add_from_rate_list(rates)
        else:
            self._rates = collections.OrderedDict()
        self._library_source_lines = collections.deque()

        if self._library_file and read_library:
            self._library_file = _find_rate_file(self._library_file)
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

    def _read_library_file(self):
        # loop through library file, read lines
        try:
            flib = open(self._library_file)
        except:
            print(f'Could not open file {self._library_file}')
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
                self._library_source_lines.popleft()
            else:
                try:
                    chapter = int(line)
                except:
                    # we can't interpret line as a chapter so use current_chapter
                    try:
                        assert(current_chapter)
                    except:
                        print(f'ERROR: malformed library file {self._library_file}, cannot identify chapter.')
                        raise
                    else:
                        chapter = current_chapter
                else:
                    self._library_source_lines.popleft()
            current_chapter = chapter

            rlines = None
            if chapter == 't':
                rlines = [self._library_source_lines.popleft() for i in range(5)]
            elif type(chapter) == int:
                rlines = [self._library_source_lines.popleft() for i in range(3)]
            if rlines:
                sio = io.StringIO('\n'.join([f'{chapter}'] +
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
        tmp_rates = [v for k, v in self._rates.items()]
        for r in sorted(tmp_rates):
            if not r.reverse:
                rstrings.append(f'{r.__repr__():30} [Q = {float(r.Q):6.2f} MeV] ({r.get_rate_id()})')
        for r in sorted(tmp_rates):
            if r.reverse:
                rstrings.append(f'{r.__repr__():30} [Q = {float(r.Q):6.2f} MeV] ({r.get_rate_id()})')

        return '\n'.join(rstrings)

    def __add__(self, other):
        """ Add two libraries to get a library containing rates from both. """
        new_rates = self._rates
        for id, r in other._rates.items():
            try:
                assert id not in new_rates
            except:
                if r != new_rates[id]:
                    print(f'ERROR: rate {r} defined differently in libraries {self._library_file} and {other._library_file}\n')
                    raise
            else:
                new_rates[id] = r
        new_library = Library(libfile=f'{self._library_file} + {other._library_file}',
                              rates=new_rates,
                              read_library=False)
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

    def get_rate(self, id):
        """ Return a rate matching the id provided. """
        try:
            return self._rates[id]
        except:
            print("ERROR: rate identifier does not match a rate in this library.")
            raise

    def diff(self, other_library):
        """Return a Library containing the rates in this library that are not
        contained in other_library"""

        diff_rates = set(self.get_rates()) - set(other_library.get_rates())
        new_library = Library(rates=diff_rates)
        return new_library

    def remove_rate(self, rate):
        """Manually remove a rate from the library by supplying the id"""

        if isinstance(rate, Rate):
            id = rate.get_rate_id()
            self._rates.pop(id)
        else:
            self._rates.pop(rate)

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

    def validate(self, other_library, forward_only=True, ostream=None):
        """perform various checks on the library, comparing to other_library,
        to ensure that we are not missing important rates.  The idea
        is that self should be a reduced library where we filtered out
        a few rates and then we want to compare to the larger
        other_library to see if we missed something important.

        ostream is the I/O stream to send output to (for instance, a
        file object or StringIO object).  If it is None, then output
        is to stdout.

        """

        current_rates = sorted(self.get_rates())

        # check the forward rates to see if any of the products are
        # not consumed by other forward rates

        passed_validation = True

        for rate in current_rates:
            if rate.reverse:
                continue
            for p in rate.products:
                found = False
                for orate in current_rates:
                    if orate == rate:
                        continue
                    if orate.reverse:
                        continue
                    if p in orate.reactants:
                        found = True
                        break
                if not found:
                    passed_validation = False
                    msg = f"validation: {p} produced in {rate} never consumed."
                    if ostream is None:
                        print(msg)
                    else:
                        ostream.write(msg + "\n")

        # now check if we are missing any rates from other_library with the exact same reactants

        for rate in current_rates:
            if forward_only and rate.reverse:
                continue

            # create a rate filter with these exact reactants
            rf = RateFilter(reactants=rate.reactants)
            all_rates_library = other_library.filter(rf)

            for other_rate in sorted(all_rates_library.get_rates()):
                # check to see if other_rate is already in current_rates
                found = True
                if other_rate not in current_rates:
                    found = False

                if not found:
                    msg = f"validation: missing {other_rate} as alternative to {rate} (Q = {other_rate.Q} MeV)."
                    if ostream is None:
                        print(msg)
                    else:
                        ostream.write(msg + "\n")

        return passed_validation

    def forward(self):
        """
        Select only the forward rates, discarding the inverse rates obtained
        by detailed balance.
        """

        only_fwd_filter = RateFilter(filter_function = lambda r: not r.reverse)
        only_fwd = self.filter(only_fwd_filter)
        return only_fwd

    def backward(self):
        """
        Select only the reverse rates, obtained by detailed balance.
        """

        only_bwd_filter = RateFilter(filter_function = lambda r: r.reverse)
        only_bwd = self.filter(only_bwd_filter)
        return only_bwd
        
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
        matches_filter_function = True
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
        if self.filter_function is not None:
            matches_filter_function = self.filter_function(r)
        return (matches_reactants and matches_products and matches_reverse and
                matches_min_reactants and matches_max_reactants and
                matches_min_products and matches_max_products and
                matches_filter_function)

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

    def __init__(self, libfile='20180319default2', rates=None, read_library=True):
        assert libfile == '20180319default2'  and rates is None and read_library, "Only the 20180319default2 default ReacLib snapshot is accepted"
        Library.__init__(self, libfile=libfile, rates=rates, read_library=read_library)
