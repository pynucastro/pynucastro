"""
Reads tabular binding energy/nucleon data file and supplies table data.
"""

from __future__ import print_function
import os

class BindingNuclide(object):
    """
    Hold the information for a single nucleus from the binding energy/nucleon table.
    """
    def __init__(self, n=None, z=None, nucbind=None):

        self.n = None
        self.z = None
        self.nucbind = None

        if n:
            self.n = int(n)
        if z:
            self.z = int(z)
        if nucbind:
            self.nucbind = float(nucbind)

    def print_contents(self):
        """
        Print Contents
        """
        print('n = {}'.format(self.n))
        print('z = {}'.format(self.z))
        print('nucbind = {}'.format(self.nucbind))

class BindingTable(object):
    """
    A simple class to manage reading and parsing the table of binding energy/nucleon.
    """

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
            fname = os.path.join(thisdir, 'AtomicMassEvaluation', 'binding_2016.txt')
            if os.path.isfile(fname):
                self.datfile = fname

        self.nuclides = []

        if self.datfile:
            self.read()

    def read(self):
        """
        Read the binding energy table
        """
        try:
            f = open(self.datfile, 'r')
        except:
            print('ERROR: data file not found!')
            exit()

        # Get rid of the header
        for _ in range(self.header_length):
            f.readline()

        # Read nuclide mass data
        for line in f:
            ls = line.strip()
            n, z, ebind = ls.split()
            nuclide = BindingNuclide(n, z, ebind)
            self.nuclides.append(nuclide)

        f.close()

    def get_nuclide(self, n=-1, z=-1):
        """
        Returns the nuclide object given n and z.
        """
        if n >= 0 and z >= 0:
            for nuc in self.nuclides:
                if nuc.n == n and nuc.z == z:
                    return nuc
        else:
            print('ERROR: invalid (n, z) supplied: ({}, {})'.format(n, z))
            exit()
        print('Nuclide not found for (n, z) = ({}, {})'.format(n, z))
