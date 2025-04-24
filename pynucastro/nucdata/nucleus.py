"""
Classes and methods to interface with files storing rate data.
"""

import re
from pathlib import Path

from pynucastro.constants import constants
from pynucastro.nucdata.elements import PeriodicTable, UnidentifiedElement
from pynucastro.nucdata.halflife_table import HalfLifeTable
from pynucastro.nucdata.mass_table import MassTable
from pynucastro.nucdata.partition_function import PartitionFunctionCollection
from pynucastro.nucdata.spin_table import SpinTable

_pynucastro_dir = Path(__file__).parents[1]
_pynucastro_rates_dir = _pynucastro_dir/'library'
_pynucastro_tabular_dir = _pynucastro_rates_dir/'tabular'

# read the various tables with nuclear properties at the module-level
_mass_table = MassTable()
_halflife_table = HalfLifeTable()
_spin_table = SpinTable(reliable=True)

# read the partition function table once and store it at the module-level
_pcollection = PartitionFunctionCollection(use_high_temperatures=True, use_set='frdm')


class UnsupportedNucleus(Exception):
    pass


class Nucleus:
    """A nucleus that participates in a reaction.

    Parameters
    ----------
    name : str
        name of the nucleus (e.g. "c12")
    dummy : bool
        a dummy nucleus is one that we can use where
        a nucleus is needed, but it is not considered
        to be part of the network

    Attributes
    ----------
    Z : float
        atomic number
    N : float
        neutron number
    A : float
        atomic mass
    short_spec_name : str
        nucleus abbreviation (e.g. "he4")
    caps_name : str
        capitalized short species name (e.g. "He4")
    el : str
        element name (e.g. "he")
    pretty : str
        LaTeX formatted version of the nucleus name
    dm : float
        mass excess (MeV)
    nucbind : float
        nuclear binding energy (MeV / nucleon)
    A_nuc : float
        nuclear mass (amu)
    mass : float
        nuclear mass (MeV)
    tau : float
        half life (s)
    spin_states : int
        the ground state spin
    partition_function : PartitionFunction
        the `PartitionFunction` object for this nucleus, which
        allows for the evaluation of the temperature-dependent
        partition function.
    dummy : bool
        is this a dummy nucleus
    nse : bool
        an NSE proton has the same properties
        as a proton but compares as being distinct
    """

    _cache = {}

    def __init__(self, name, dummy=False):
        name = name.lower()
        self.raw = name

        self.dummy = dummy
        self.nse = False

        # element symbol and atomic weight
        if name == "p":
            self.el = "h"
            self.A = 1
            self.short_spec_name = "h1"
            self.caps_name = "p"
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
            self.el = "he"
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
            self.caps_name = "n"
        elif name == "p_nse":
            # this is a proton with a different name
            # it is meant to be used in iron-group rates
            # in NSE
            self.el = "h"
            self.A = 1
            self.Z = 1
            self.N = 0
            self.short_spec_name = "p_nse"
            self.spec_name = "proton-nse"
            self.pretty = r"\mathrm{p}_\mathrm{NSE}"
            self.caps_name = "p_NSE"
            self.nse = True
        elif name.strip() in ("al-6", "al*6"):
            raise UnsupportedNucleus()
        else:
            e = re.match(r"([a-zA-Z]*)(\d*)", name)
            self.el = e.group(1).title()  # chemical symbol
            assert self.el
            self.A = int(e.group(2))
            assert self.A >= 0
            self.short_spec_name = name
            self.caps_name = name.capitalize()

        # use lowercase element abbreviation regardless the case of the input
        self.el = self.el.lower()

        # atomic number comes from periodic table
        if name not in ["n", "p_nse"]:
            i = PeriodicTable.lookup_abbreviation(self.el)
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

        # set the number of spin states
        try:
            self.spin_states = _spin_table.get_spin_states(a=self.A, z=self.Z)
        except NotImplementedError:
            self.spin_states = None

        # set a partition function object to every nucleus
        try:
            self.partition_function = _pcollection.get_partition_function(self.short_spec_name)
        except ValueError:
            self.partition_function = None

        # nuclear mass
        try:
            # Note: for Nubase 2020, we need to use the CODATA 18 constants
            mass_H = _mass_table.get_mass_diff(a=1, z=1) + constants.m_u_MeV_C18
            self.dm = _mass_table.get_mass_diff(a=self.A, z=self.Z)
            self.A_nuc = float(self.A) + self.dm / constants.m_u_MeV_C18
            self.mass = self.A * constants.m_u_MeV_C18 + self.dm
            B = (self.Z * mass_H + self.N * constants.m_n_MeV_C18) - self.mass
            self.nucbind = B / self.A

        except NotImplementedError:
            self.dm = None
            self.A_nuc = None
            self.mass = None
            self.nucbind = None

        # halflife
        try:
            self.tau = _halflife_table.get_halflife(a=self.A, z=self.Z)
        except NotImplementedError:
            self.tau = None

    @classmethod
    def from_cache(cls, name, dummy=False):
        """Check if we've already created this nucleus, and if so,
        return a reference to it from the cache.

        Parameters
        ----------
        name : str
            name of the nucleus (e.g. "c12")
        dummy : bool
            a dummy nucleus is one that we can use where
            a nucleus is needed, but it is not considered
            to be part of the network

        Returns
        -------
        Nucleus
        """

        key = (name.lower(), dummy)
        if key not in cls._cache:
            cls._cache[key] = cls(name, dummy)
        return cls._cache[key]

    @classmethod
    def from_Z_A(cls, Z, A, dummy=False):
        """Create a nucleus given Z and A

        Parameters
        ----------
        Z : int
            atomic number
        A : int
            atomic weight
        dummy : bool
            a dummy nucleus is one that we can use where
            a nucleus is needed, but it is not considered
            to be part of the network

        Returns
        -------
        Nucleus
        """

        # checks if Z and A are valid inputs
        if not (isinstance(Z, int) and isinstance(A, int)):
            raise TypeError("Nucleus Z and A must be integers")
        if not (Z >= 0 and A >= 0):
            raise ValueError("Nucleus Z and A must be non-negative")
        if Z > A:
            raise ValueError("Nucleus Z can't be bigger than A")

        # checks if neutron
        if (Z, A) == (0, 1):
            return cls.from_cache("n", dummy)

        # otherwise, finds element Z on the periodic table
        i = PeriodicTable.lookup_Z(Z)
        if i is None:
            raise UnidentifiedElement(f"Element {Z} could not be found")

        name = i.abbreviation + str(A)
        return cls.from_cache(name, dummy)

    def summary(self):
        """print a summary of the nuclear properties"""

        print(f"{self.caps_name} / {self.spec_name}:")
        print(f"    A: {self.A}")
        print(f"    N: {self.N}")
        print(f"    Z: {self.Z}")

        print("")
        print(f"    mass: {self.mass:.5f} MeV")
        print(f"    mass excess: {self.dm:.5f} MeV")
        print(f"    binding energy: {self.nucbind:.5f} MeV")

        print("")
        if self.tau == "stable":
            print(f"    half-life: {self.tau}")
        else:
            print(f"    half-life: {self.tau} s")

        print("")
        if self.partition_function:
            print("    partition function: available")
        else:
            print("    partition function: not available")

        if self.spin_states:
            print(f"    spin states: {self.spin_states}")
        else:
            print("    spin states: not available")

        print("")
        print(f"    dummy: {self.dummy}")
        print(f"    nse: {self.nse}")

    def __repr__(self):
        if self.raw not in ("p", "p_nse", "d", "t", "n"):
            return self.raw.capitalize()
        return self.raw

    def __hash__(self):
        return hash((self.Z, self.A))

    def c(self):
        """Return the capitalized-style name"""
        return self.caps_name

    def cindex(self):
        """Return the name for C++ indexing"""
        return self.short_spec_name.capitalize()

    def __eq__(self, other):
        if isinstance(other, Nucleus):
            return (self.el == other.el and
                    self.Z == other.Z and self.A == other.A and
                    self.nse == other.nse)
        if isinstance(other, tuple):
            return (self.Z, self.A) == other
        return NotImplemented

    def __lt__(self, other):
        if not self.Z == other.Z:
            return self.Z < other.Z
        return self.A < other.A

    def __add__(self, other):
        Z = self.Z + other.Z
        A = self.A + other.A
        dummy = self.dummy and other.dummy
        return Nucleus.from_Z_A(Z, A, dummy)

    def __sub__(self, other):
        Z = self.Z - other.Z
        A = self.A - other.A
        dummy = self.dummy and other.dummy
        return Nucleus.from_Z_A(Z, A, dummy)

    @classmethod
    def cast(cls, obj):
        """Create a Nucleus from a string

        Parameters
        ----------
        obj : str or Nucleus
            the object to cast.  If it is a `Nucleus`, then we simply return
            it.

        Returns
        -------
        Nucleus
        """
        if isinstance(obj, cls):
            return obj
        if isinstance(obj, str):
            return cls.from_cache(obj)
        raise TypeError("Invalid type passed to Nucleus.cast() (expected str or Nucleus)")

    @classmethod
    def cast_list(cls, lst, *, allow_None=False, allow_single=False):
        """Convert a list of objects into a list of Nucleus objects

        Parameters
        ----------
        lst : list
            a list of `str` or `Nucleus`
        allow_None : bool
            allow lst = None and simply return None
        allow_single : bool
            allow lst to be a single `str` or `Nucleus` instead
            of a `list`

        Returns
        -------
        list
        """

        if allow_None and lst is None:
            return lst
        if isinstance(lst, (str, cls)):
            if allow_single:
                return [cls.cast(lst)]
            raise ValueError("Single object passed to Nucleus.cast_list() instead of list")
        return [cls.cast(obj) for obj in lst]


def get_nuclei_in_range(zmin, zmax, amin, amax):
    """Given a range of Z = [zmin, zmax], and A = [amin, amax],
    return a list of Nucleus objects for all nuclei in this range

    Parameters
    ----------
    zmin : int
        minimum atomic number
    zmax : int
        maximum atomic number
    amin : int
        minimum atomic weight
    amax : int
        maximum atomic weight

    Returns
    -------
    list
    """

    nuc_list = []
    assert zmax >= zmin, "zmax must be >= zmin"
    assert amax >= amin, "amax must be >= amin"

    for z in range(zmin, zmax+1):
        element = PeriodicTable.lookup_Z(z)
        for a in range(amin, amax+1):
            name = f"{element.abbreviation}{a}"
            nuc_list.append(Nucleus(name))

    return nuc_list


def get_all_nuclei():
    """Return a list with every Nucleus that has a known mass

    Returns
    -------
    list
    """

    nuc_list = []

    for (A, Z) in _mass_table.mass_diff:
        if Z == 0 and A == 1:
            nuc = "n"
        else:
            el = PeriodicTable.lookup_Z(Z)
            nuc = f"{el.abbreviation}{A}"
        nuc_list.append(Nucleus(nuc))

    return nuc_list
