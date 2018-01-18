# Common Imports
from __future__ import print_function

class BindingNuclide(object):
    """Hold the information for a single nucleus from the binding energy/nucleon table."""
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
