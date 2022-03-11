"""
Classes and methods to interface with files storing rate data.
"""

import os
import re
import io
import collections

import numpy as np
import matplotlib.pyplot as plt
import numba

try:
    from numba.experimental import jitclass
except ImportError:
    from numba import jitclass

from pynucastro.nucdata import UnidentifiedElement, PeriodicTable

_pynucastro_dir = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
_pynucastro_rates_dir = os.path.join(_pynucastro_dir, 'library')
_pynucastro_tabular_dir = os.path.join(_pynucastro_rates_dir, 'tabular')


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
    """ precompute temperature factors for speed """

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
    """

    def __init__(self, a, labelprops=None):
        """here a is iterable (e.g., list or numpy array), storing the
           coefficients, a0, ..., a6

        """
        self.a = a
        self.labelprops = labelprops
        self._update_label_properties()
        self.label = None
        self.resonant = None
        self.weak = None
        self.reverse = None

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
    pretty printing string

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

        # use lowercase element abbreviation regardless the case of the input
        self.el = self.el.lower()

        # set a partition function object to every nucleus
        self._partition_function = None

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

    def set_partition_function(self, p_collection, set_data='frdm', use_high_temperatures=True):
        """
        This function associates to every nucleus a PartitionFunction object.
        """
        assert type(p_collection) == PartitionFunctionCollection

        p_collection.set_data_selector(set_data)
        p_collection.use_high_temperatures(use_high_temperatures)
        self._partition_function = p_collection.get_partition_function(self)

    def get_partition_function(self):
        return self._partition_function

    def __repr__(self):
        return self.raw

    def __hash__(self):
        return hash((self.Z, self.A))

    def c(self):
        return self.caps_name

    def __eq__(self, other):
        if isinstance(other, Nucleus):
            return self.el == other.el and \
               self.Z == other.Z and self.A == other.A
        elif isinstance(other, tuple):
            return (self.Z, self.A) == other
        else:
            return NotImplemented

    def __lt__(self, other):
        if not self.Z == other.Z:
            return self.Z < other.Z
        else:
            return self.A < other.A

