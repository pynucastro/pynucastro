import collections
import io
import re
from os import walk
from pathlib import Path

from pynucastro.nucdata import Nucleus, UnsupportedNucleus
from pynucastro.rates.derived_rate import DerivedRate
from pynucastro.rates.files import (RateFileError, _find_rate_file,
                                    get_rates_dir)
from pynucastro.rates.known_duplicates import (find_duplicate_rates,
                                               is_allowed_dupe)
from pynucastro.rates.rate import Rate
from pynucastro.rates.reaclib_rate import ReacLibRate
from pynucastro.rates.tabular_rate import TabularRate


def list_known_rates():
    """Print a list of all of the rates found in the library """

    lib_path = Path(__file__).parents[1]/"library"

    for _, _, filenames in walk(lib_path):
        for f in filenames:
            # skip over files that are not rate files
            if Path(f).suffix in (".md", ".dat", ".py", "ipynb"):
                continue
            try:
                lib = Library(f)
            except (RateFileError, UnsupportedNucleus):
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
            if nuc.lower() == "pp":
                reactants += [Nucleus("p"), Nucleus("p")]
                continue
            if nuc.lower() == "aa":
                reactants += [Nucleus("he4"), Nucleus("he4")]
                continue
            print(f"couldn't deal with {nuc}")
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
            if nuc.lower() == "pp":
                products += [Nucleus("p"), Nucleus("p")]
                continue
            if nuc.lower() == "aa":
                products += [Nucleus("he4"), Nucleus("he4")]
                continue
            print(f"couldn't deal with {nuc}")
            raise

    return reactants, products


def capitalize_rid(rid, delimiter):
    # Used to capitalize rid or fname given the delimiter
    # delimiter is usually either "_" or " "

    rid_nucs = rid.split(delimiter)
    rid_mod = []
    do_capitalization = True
    for n in rid_nucs:
        if n in ("weak", "approx", "derived"):
            do_capitalization = False
        if do_capitalization:
            if n not in ("n", "p"):
                n = n.capitalize()
        rid_mod.append(n)

    rid_mod = delimiter.join(rid_mod)
    return rid_mod


class Library:
    """A Library is a container storing multiple rates that allows for
    filtering rates based on rules, managing duplicate rates, and
    selecting subsets of rates based on properties.  At its heart is a
    ``dict`` of rates keyed by the rate id.

    A library may contain rates from a single source, or be created by
    adding or subtracting existing Library objects.

    Parameters
    ----------
    libfile : str
        a file containing a sequence of rates in a format that we
        understand (for example a ReacLib database)
    rates : list, dict, Rate
        a single :py:class:`Rate <pynucastro.rates.rate.Rate>` or an
        iterable of `Rate` objects.  If it is a dictionary, then it
        should be keyed by the rate id.

    """

    def __init__(self, libfile=None, rates=None):
        self._rates = {}

        if rates:
            if isinstance(rates, Rate):
                rates = [rates]
            if isinstance(rates, dict):
                self._rates = rates
            elif isinstance(rates, (list, set)):
                self.add_rates(rates)
            else:
                raise TypeError("rates in Library constructor must be a Rate object, list of Rate objects, or dictionary of Rate objects keyed by Rate.id")

        if libfile:
            library_file = _find_rate_file(libfile)
            self._read_library_file(library_file)

    def get_rates(self):
        """Return a list of the rates in this library.

        Returns
        -------
        list

        """
        return list(self._rates.values())

    def get_rate(self, rid):
        """Return a rate matching the id provided.

        Returns
        -------
        Rate

        """

        try:
            rid_mod = capitalize_rid(rid, " ")
            return self._rates[rid_mod]
        except KeyError:
            pass

        # fallback to the rate fname
        try:
            rid_mod = capitalize_rid(rid, "_")
            return [q for q in self.get_rates() if q.fname == rid_mod][0]
        except IndexError:
            raise LookupError(f"rate identifier {rid!r} does not match a rate in this library.") from None

    @property
    def num_rates(self):
        """Get the total number of rates in the Library

        Returns
        -------
        int

        """
        return len(self.get_rates())

    def add_rate(self, rate):
        """Manually add a rate to the library.

        Parameters
        ----------
        rate : Rate
            The rate to add

        """

        if not isinstance(rate, Rate):
            raise TypeError(f"invalid Rate object {rate}")
        rid = rate.id

        if rid in self._rates:
            raise ValueError(f"supplied a Rate object already in the Library: {rid}")
        self._rates[rid] = rate

    def add_rates(self, ratelist):
        """Add multiple rates to the library

        Parameters
        ----------
        ratelist : list of Rate
            the list of rates to add to the library.

        """

        for rate in ratelist:
            self.add_rate(rate)

    def get_rate_by_name(self, name):
        """Given a string representing a rate in the form 'A(x,y)B'
        (or a list of strings for multiple rates) return the Rate
        objects that match from the Library.  If there are multiple
        inputs, then a list of Rate objects is returned.

        Returns
        -------
        rates : list, Rate
            A single rate or a list of rates

        """

        rate_name_list = name
        if isinstance(name, str):
            rate_name_list = [name]

        rates_out = []

        for rname in rate_name_list:
            reactants, products = _rate_name_to_nuc(rname)

            rf = RateFilter(reactants=reactants, products=products)
            _lib = self.filter(rf)
            if _lib is None:
                print(f"rate {rname} not found")
                continue
            rates_out += _lib.get_rates()

        if len(rates_out) == 0:
            return None
        if len(rates_out) == 1:
            return rates_out[0]
        return rates_out

    def remove_rate(self, rate):
        """Manually remove a rate from the library by supplying the
        short name "A(x,y)B", a Rate object, or the rate id

        Parameters
        ----------
        rate : str, Rate
            The rate to remove from the library.

        """

        if isinstance(rate, Rate):
            rid = rate.id
            self._rates.pop(rid)
        elif isinstance(rate, str):
            rid = self.get_rate_by_name(rate).id
            self._rates.pop(rid)
        else:
            # we assume that a rate id as provided
            self._rates.pop(rate)

    def get_nuclei(self):
        """Get the list of unique nuclei in the library

        Returns
        -------
        set

        """

        return {nuc for r in self.get_rates() for nuc in r.reactants + r.products}

    def heaviest(self):
        """Return the heaviest nuclide in this library.

        Returns
        -------
        Nucleus

        """

        nuc = None
        for r in self.get_rates():
            rnuc = r.heaviest()
            if nuc:
                if rnuc.A > nuc.A or (rnuc.A == nuc.A and rnuc.Z < nuc.Z):
                    nuc = rnuc
            else:
                nuc = rnuc
        return nuc

    def lightest(self):
        """Return the lightest nuclide in this library.

        Returns
        -------
        Nucleus

        """

        nuc = None
        for r in self.get_rates():
            rnuc = r.lightest()
            if nuc:
                if rnuc.A < nuc.A or (rnuc.A == nuc.A and rnuc.Z > nuc.Z):
                    nuc = rnuc
            else:
                nuc = rnuc
        return nuc

    def _read_library_file(self, library_file):
        # loop through library file, read lines

        library_source_lines = collections.deque()

        with library_file.open("r") as flib:
            for line in flib:
                ls = line.rstrip('\n')
                if ls.strip():
                    library_source_lines.append(ls)

        # identify distinct rates from library lines
        current_chapter = None
        while True:
            if len(library_source_lines) == 0:
                break

            # Check to see if there is a chapter ID, if not then use current_chapter
            # (for Reaclib v1 formatted library files)
            line = library_source_lines[0].strip()
            chapter = None
            if line in ('t', 'T'):
                chapter = 't'
                library_source_lines.popleft()
            else:
                try:
                    chapter = int(line)
                except (TypeError, ValueError):
                    # we can't interpret line as a chapter so use current_chapter
                    assert current_chapter, f'malformed library file {library_file}, cannot identify chapter.'
                    chapter = current_chapter
                else:
                    library_source_lines.popleft()
            current_chapter = chapter

            rlines = None
            rate_type = None
            if chapter == 't':
                rlines = [library_source_lines.popleft() for i in range(5)]
                rate_type = "tabular"
            elif isinstance(chapter, int):
                rlines = [library_source_lines.popleft() for i in range(3)]
                rate_type = "reaclib"
            if rlines:
                sio = io.StringIO('\n'.join([f'{chapter}'] +
                                            rlines))
                #print(sio.getvalue())
                try:
                    if rate_type == "reaclib":
                        r = ReacLibRate(rfile=sio)
                    elif rate_type == "tabular":
                        r = TabularRate(rfile=sio)
                    else:
                        raise NotImplementedError("rate not implemented")
                except UnsupportedNucleus:
                    pass
                else:
                    rid = r.id
                    if rid in self._rates:
                        self._rates[rid] = self._rates[rid] + r
                    else:
                        self._rates[rid] = r

    def write_to_file(self, filename, *, prepend_rates_dir=False):
        """Write the library out to a file of the given name in
        Reaclib format.

        Parameters
        ----------
        filename : str
            The filename to use for the library
        prepend_rates_dir : bool
            If ``True``, then output to the pynucastro rate file
            directory.

        """

        if prepend_rates_dir:
            filename = get_rates_dir()/filename

        with filename.open("w") as f:
            for rate in self.get_rates():
                rate.write_to_file(f)

    def __repr__(self):
        """Return a string containing the rates IDs in this library."""
        rstrings = []
        tmp_rates = [v for k, v in self._rates.items()]
        for r in sorted(tmp_rates):
            if r.Q is not None and r.Q >= 0:
                rstrings.append(f'{r.__repr__():30} [Q = {float(r.Q):6.2f} MeV] ({r.id})')
        for r in sorted(tmp_rates):
            if r.Q is not None and r.Q < 0:
                rstrings.append(f'{r.__repr__():30} [Q = {float(r.Q):6.2f} MeV] ({r.id})')

        for r in sorted(tmp_rates):
            if r.Q is None:
                rstrings.append(f'{r.__repr__():30} ({r.id})')

        return '\n'.join(rstrings)

    def __add__(self, other):
        """Add two libraries to get a library containing rates from
        both.

        """
        new_rates = self._rates
        for rid, r in other._rates.items():
            if rid in new_rates:
                if r != new_rates[rid]:
                    raise ValueError(f'rate {r} defined differently in libraries')
            else:
                new_rates[rid] = r
        new_library = Library(rates=new_rates)
        return new_library

    def __sub__(self, other):
        """Return a Library containing the rates in this library that
        are not contained in other_library

        """

        diff_rates = set(self.get_rates()) - set(other.get_rates())
        new_library = Library(rates=diff_rates)
        return new_library

    def get_rate_by_nuclei(self, reactants, products):
        """Given a list of reactants and products, return any matching
        rates

        Parameters
        ----------
        reactants : list of Nucleus or str
            the list of nuclei that serve as reactants.
        products : list of Nucleus or str
            the list of nuclei that serve as products.

        Returns
        -------
        Rate, list(Rate)
            a list of Rate object or a single Rate (if there is only one)

        """

        reactants = sorted(Nucleus.cast_list(reactants))
        products = sorted(Nucleus.cast_list(products))
        _tmp = [r for r in self.get_rates() if
                sorted(r.reactants) == reactants and
                sorted(r.products) == products]

        if not _tmp:
            return None
        if len(_tmp) == 1:
            return _tmp[0]
        return _tmp

    def find_duplicate_links(self):
        """Find instances of multiple rates having the same reactants
        and products.  These may not be the same Rate object (e.g.,
        one could be tabular the other a simple decay), but they will
        present themselves in the network as the same link.

        Returns
        -------
        duplicate_rates : list
            a list where each entry is a list of all the rates
            that share the same link.

        """

        duplicates = find_duplicate_rates(self.get_rates())

        # there are some allowed duplicates for special cases.  We
        # will now check for those
        dupe_to_remove = []
        for dupe in duplicates:
            if is_allowed_dupe(dupe):
                dupe_to_remove.append(dupe)

        for dupe in dupe_to_remove:
            duplicates.remove(dupe)

        return duplicates

    def linking_nuclei(self, nuclist, *, with_reverse=True,
                       print_warning=True):
        """Return a library containing the rates linking the list of
        nuclei passed in.

        Parameters
        ----------
        nuclist : list of str or Nucleus
            the nuclei to link (either the string names or the Nucleus objects)
        with_reverse : bool
            do we include reverse rates?
        print_warning : bool
            if ``True``, then print a warning if one of the input
            nuclei is not linked.

        Returns
        -------
        Library

        """

        nucleus_set = set(Nucleus.cast_list(nuclist))

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
        if print_warning:
            lib_nuclei = new_lib.get_nuclei()
            for nuc in nucleus_set:
                if nuc not in lib_nuclei:
                    print(f"warning: {nuc} was not able to be linked")

        return new_lib

    def filter(self, filter_spec):
        """Filter the rates in the library based on a set of rules.

        Parameters
        ----------
        filter_spec : RateFilter, list
            a filter (or list of filters) to apply to the library
            to define a subset of rates.

        Returns
        -------
        Library

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
        """Select only the forward rates, discarding the inverse rates
        obtained by detailed balance.  Note: "forward" here means
        that it is not a reverse rate derived from detailed balance,
        and does not necessarily mean Q > 0.

        Returns
        -------
        Library

        """

        only_fwd_filter = RateFilter(reverse=False)
        only_fwd = self.filter(only_fwd_filter)
        return only_fwd

    def backward(self):
        """Select only the reverse rates, obtained by detailed
        balance.  Note: "reverse" here means that it was derived
        by detailed balance, and not that Q < 0.

        Returns
        -------
        Library

        """

        only_bwd_filter = RateFilter(reverse=True)
        only_bwd = self.filter(only_bwd_filter)
        return only_bwd

    def forward_for_detailed_balance(self):
        """Loop over the forward rates (as filtered by
        :py:meth:`.forward`) and return those that can be used to
        derive reverse rates via detailed balance.  This means that
        they cannot be tabular or weak rates.

        Returns
        -------
        Library

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
        """Loop over all of the forward rates that can be used to
        derive inverse rates (as returned by
        :py:meth:`.forward_for_detailed_balance`) and derive the
        inverses, potentially taking into account the partition
        function and recomputing Q.

        Parameters
        ----------
        compute_Q : bool
            do we recompute the Q value based on the masses?
        use_pf : bool
            do we use the temperature-dependent partition function?

        Returns
        -------
        Library

        """

        derived_rates = []
        onlyfwd = self.forward_for_detailed_balance()

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
            self.reactants = Nucleus.cast_list(reactants, allow_single=True)
        if products:
            self.products = Nucleus.cast_list(products, allow_single=True)

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
    """Create a :py:class:`Library` containing all of the rates in the
    latest stored version of the ReacLib library.

    """

    def __init__(self):
        libfile = 'reaclib_default2_20250330'
        Library.__init__(self, libfile=libfile)


class TabularLibrary(Library):
    """Create a :py:class:`Library` containing all of the tabular
    rates we know (excluding duplications) across multiple sources.

    Parameters
    ----------
    ordering : list of str
        The list of sources of the rates from lowest to highest
        precedence.  We will read from the first source, and then for
        any later sources, for any duplicate rates, we will replace
        the existing rate with the version from the higher-priority
        library.  The default ordering is ``["ffn", "langanke",
        "suzuki"]``

    """

    lib_path = Path(__file__).parents[1]/"library/tabular"

    def __init__(self, ordering=None):
        # find all of the tabular rates that pynucastro knows about
        # we'll assume that these are of the form *betadecay.dat or
        # *electroncapture.dat

        if ordering is None:
            ordering = ["ffn", "langanke", "suzuki"]

        trates = []

        for source in ordering:
            source_dir = self.lib_path / Path(source)
            for _, _, filenames in sorted(walk(source_dir)):
                for f in sorted(filenames):
                    if f.endswith("electroncapture.dat") or f.endswith("betadecay.dat"):
                        r = TabularRate(rfile=source_dir / f)
                        if r in trates:
                            # we are looping over the various libraries in order
                            # from lowest precedence to highest.  So if the rate
                            # exists, then delete it and add this one.  Since
                            # matching only looks at reactants and products, we
                            # can use the new r for both deleting a adding
                            trates.remove(r)
                        trates.append(r)

        Library.__init__(self, rates=trates)


class SuzukiLibrary(TabularLibrary):
    """Create a :py:class:`Library` containing all of the tabular
    rates inside the "suzuki" subdirectory.

    """

    def __init__(self):
        super().__init__(ordering=["suzuki"])


class LangankeLibrary(TabularLibrary):
    """Create a :py:class:`Library` containing all of the tabular
    rates inside the "langanke" subdirectory.

    """

    def __init__(self):
        super().__init__(ordering=["langanke"])


class FFNLibrary(TabularLibrary):
    """Create a :py:class:`Library` containing all of the tabular
    rates inside the "ffn" subdirectory.

    """

    def __init__(self):
        super().__init__(ordering=["ffn"])
