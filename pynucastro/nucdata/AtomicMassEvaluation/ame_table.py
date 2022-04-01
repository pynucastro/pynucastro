"""
Reads AME tabular mass data file and supplies table data.
"""

# Common Imports
import os

from ame_nuclide import AMENuclide


def str_head_pop(s, n):
    """
    Pop n characters from string s from the beginning.

    Returns a tuple consisting of the new s and the popped characters.
    """
    if n >= len(s):
        return ('', s)
    s = list(s)
    return (''.join(s[n:]), ''.join(s[:n]))


def cleanup(s):
    """
    Given string s, removes "#" and "*".

    # is a standin for a decimal point for estimated values.

    * represents a non-computable quantity.
    """
    s = s.replace('#', '.')
    s = s.replace('*', '')
    s = s.replace(' ', '')
    return s


class AMETable(object):
    """A simple class to manage reading and parsing the AME Table data
    files (2012, 2016)."""

    header_length = 39

    def __init__(self, datfile=None):
        """
        Initialize
        """
        self.datfile = None
        if datfile:
            self.datfile = datfile
        else:
            thisdir = os.path.dirname(os.path.realpath(__file__))
            fname = os.path.join(thisdir, 'mass.mas12')
            if os.path.isfile(fname):
                self.datfile = fname

        self.nuclides = []

        if self.datfile:
            self.read()

    def read(self):
        """
        Read the AME data file 'datfile'
        """
        try:
            f = open(self.datfile, 'r')
        except IOError:
            print('ERROR: data file not found!')
            exit()

        # Get rid of the header
        for _ in range(self.header_length):
            f.readline()

        # Read nuclide mass data
        for line in f:
            ls = line.strip()
            if ls:
                # Pull data out of the line
                line, _ = str_head_pop(line, 4)
                line, n = str_head_pop(line, 5)
                line, z = str_head_pop(line, 5)
                line, a = str_head_pop(line, 5)
                line, _ = str_head_pop(line, 1)
                line, element = str_head_pop(line, 3)
                line, origin = str_head_pop(line, 4)
                line, _ = str_head_pop(line, 1)
                line, mexcess = str_head_pop(line, 13)
                mexcess = cleanup(mexcess)
                line, d_mexcess = str_head_pop(line, 11)
                d_mexcess = cleanup(d_mexcess)
                line, nucbind = str_head_pop(line, 11)
                nucbind = cleanup(nucbind)
                line, d_nucbind = str_head_pop(line, 9)
                d_nucbind = cleanup(d_nucbind)
                line, _ = str_head_pop(line, 1)
                line, decay_type = str_head_pop(line, 2)
                line, ebeta = str_head_pop(line, 11)
                ebeta = cleanup(ebeta)
                line, d_ebeta = str_head_pop(line, 9)
                d_ebeta = cleanup(d_ebeta)
                line, _ = str_head_pop(line, 1)
                line, mass = str_head_pop(line, 16)
                mass = cleanup(mass)
                line, d_mass = str_head_pop(line, 11)
                d_mass = cleanup(d_mass)

                # Store nuclide data
                nuclide = AMENuclide(n, z, a, element.strip(), origin.strip(),
                                     mexcess, d_mexcess, nucbind,
                                     d_nucbind, decay_type.strip(), ebeta,
                                     d_ebeta, mass, d_mass)
                self.nuclides.append(nuclide)
        f.close()

    def abbrev_get_nuclide(self, isostring):
        """
        Returns the nuclide object given an identifying string
        isostring which consists of:

        [A]-[Abbreviation]

        E.g. "4-He", "7-li", etc.
        """
        try:
            a, abb = tuple(isostring.split('-'))
        except ValueError:
            print('ERROR: Supply an isotope string in the form "4-He"')
            raise
        a = int(a)
        abb = abb.strip().lower()
        for nuc in self.nuclides:
            if nuc.a == a and nuc.element.lower() == abb:
                return nuc
        print('ERROR: Could not find a nuclide by the specification {}'.format(isostring))

    def get_nuclide(self, n=-1, z=-1, a=-1):
        """
        Returns the nuclide object given at least 2 of n, z, a.

        Complains if A-Z-N != 0.
        """
        if n >= 0 and z >= 0 and a >= 0:
            if a - z - n != 0:
                print('ERROR: Z+N != A')
                exit()
        else:
            if n == -1:
                n = a - z
            elif z == -1:
                z = a - n
            else:
                a = n + z

        for nuc in self.nuclides:
            if nuc.n == n and nuc.z == z:
                return nuc

        print('Nuclide not found for (N, Z, A) = ({}, {}, {})!'.format(n, z, a))
