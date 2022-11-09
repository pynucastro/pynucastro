"""
Classes and methods to interface with files storing rate data.
"""

import os
import re

from scipy.constants import physical_constants

from pynucastro.nucdata.binding_table import BindingTable
from pynucastro.nucdata.elements import PeriodicTable
from pynucastro.nucdata.mass_table import MassTable
from pynucastro.nucdata.partition_function import PartitionFunctionCollection
from pynucastro.nucdata.spin_table import SpinTable

_pynucastro_dir = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
_pynucastro_rates_dir = os.path.join(_pynucastro_dir, 'library')
_pynucastro_tabular_dir = os.path.join(_pynucastro_rates_dir, 'tabular')

#set the atomic mass unit constant in MeV
m_u, _, _ = physical_constants['atomic mass constant energy equivalent in MeV']

#read the mass excess table once and store it at the module-level
_mass_table = MassTable()

#read the spin table once and store it at the module-level
_spin_table = SpinTable(reliable=True)

# read the binding energy table once and store it at the module-level
_binding_table = BindingTable()

# read the partition function table once and store it at the module-level
_pcollection = PartitionFunctionCollection(use_high_temperatures=True, use_set='frdm')


class UnsupportedNucleus(Exception):
    pass


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
    :var A_nuc:           Nuclear Mass in amu

    """
    _cache = {}

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
        if name != "n":
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

        try:
            self.nucbind = _binding_table.get_binding_energy(n=self.N, z=self.Z)
        except NotImplementedError:
            # the binding energy table doesn't know about this nucleus
            self.nucbind = None

        # Now we will define the Nuclear Mass,
        try:
            self.A_nuc = float(self.A) + _mass_table.get_mass_diff(a=self.A, z=self.Z) / m_u
        except NotImplementedError:
            self.A_nuc = None

    @classmethod
    def from_cache(cls, name, dummy=False):
        key = (name.lower(), dummy)
        if key not in cls._cache:
            cls._cache[key] = Nucleus(name, dummy)
        return cls._cache[key]

    def __repr__(self):
        return self.raw

    def __hash__(self):
        return hash((self.Z, self.A))

    def c(self):
        """return the capitalized-style name"""
        return self.caps_name

    def cindex(self):
        """return the name for C++ indexing"""
        return self.short_spec_name.capitalize()

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
