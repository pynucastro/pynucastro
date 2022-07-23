"""
Reads tabular binding energy/nucleon data file and supplies table data.
"""

# Common Imports
import os

from pynucastro.nucdata.binding_nuclide import BindingNuclide


class BindingTable:
    """A simple class to manage reading and parsing the table of binding energy/nucleon."""

    header_length = 2

    def __init__(self, datfile=None):
        """
        Initialize
        """
        self.datfile = None
        if datfile:
            self.datfile = datfile
        else:
            thisdir = os.path.dirname(os.path.realpath(__file__))
            fname = os.path.join(thisdir, 'binding_2016.txt')
            if os.path.isfile(fname):
                self.datfile = fname

        self.nuclides = {}

        if self.datfile:
            self.read()

    def read(self):
        """
        Read the binding energy table
        """
        try:
            f = open(self.datfile, 'r')
        except IOError:
            print('ERROR: data file not found!')
            raise

        # Get rid of the header
        for _ in range(self.header_length):
            f.readline()

        # Read nuclide mass data
        for line in f:
            ls = line.strip()
            n, z, ebind = ls.split()
            nuclide = BindingNuclide(n, z, ebind)
            self.nuclides[f"{n}_{z}"] = nuclide

        f.close()

    def get_nuclide(self, n=-1, z=-1):
        """
        Returns the nuclide object given n and z.
        """
        try:
            return self.nuclides[f"{n}_{z}"]
        except KeyError:
            raise NotImplementedError(f"nuclear data for Z={z} and N={n} not available") from None
