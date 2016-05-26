"""
Reads the AME 2012 mass data file and supplies table data.
"""

from __future__ import print_function
import numpy as np

def str_head_pop(s,n):
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
    Given string s, removes # and *.
    """
    s = s.replace('#','')
    s = s.replace('*','')
    s = s.replace(' ','')
    return s

class Nuclide(object):
    def __init__(self, n=None, z=None, a=None, element=None, origin=None,
                 mexcess=None, d_mexcess=None, nucbind=None,
                 d_nucbind=None, decay_type=None, ebeta=None,
                 d_ebeta=None, mass=None, d_mass=None):

        self.n=None
        self.z=None
        self.a=None
        self.element=None
        self.origin=None
        self.mexcess=None
        self.d_mexcess=None
        self.nucbind=None
        self.d_nucbind=None
        self.decay_type=None
        self.ebeta=None
        self.d_ebeta=None
        self.mass=None
        self.d_mass=None
        
        if n:
            self.n = int(n)
        if z:
            self.z = int(z)
        if a:
            self.a = int(a)
        self.element = element
        self.origin = origin
        if mexcess:
            self.mexcess = float(mexcess)
        if d_mexcess:
            self.d_mexcess = float(d_mexcess)
        if nucbind:
            self.nucbind = float(nucbind)
        if d_nucbind:
            self.d_nucbind = float(d_nucbind)
        self.decay_type = decay_type
        if ebeta:
            self.ebeta = float(ebeta)
        if d_ebeta:
            self.d_ebeta = float(d_ebeta)
        if mass:
            self.mass = float(mass)
        if d_mass:
            self.d_mass = float(d_mass)

        self.convert_MeV()
        self.convert_amu()

    def print_contents(self):
        """
        Print Contents
        """
        print('n = {}'.format(self.n))
        print('z = {}'.format(self.z))
        print('a = {}'.format(self.a))
        print('element = {}'.format(self.element))
        print('origin = {}'.format(self.origin))
        print('mexcess = {}'.format(self.mexcess))
        print('d_mexcess = {}'.format(self.d_mexcess))
        print('nucbind = {}'.format(self.nucbind))
        print('d_nucbind = {}'.format(self.d_nucbind))
        print('decay_type = {}'.format(self.decay_type))
        print('ebeta = {}'.format(self.ebeta))
        print('d_ebeta = {}'.format(self.d_ebeta))
        print('mass = {}'.format(self.mass))
        print('d_mass = {}'.format(self.d_mass))
        
    def convert_MeV(self):
        """
        Convert keV to MeV
        """
        if self.mexcess:
            self.mexcess = self.mexcess/1.0e3
        if self.d_mexcess:
            self.d_mexcess = self.d_mexcess/1.0e3
        if self.nucbind:
            self.nucbind = self.nucbind/1.0e3
        if self.d_nucbind:
            self.d_nucbind = self.d_nucbind/1.0e3
        if self.ebeta:
            self.ebeta = self.ebeta/1.0e3
        if self.d_ebeta:
            self.d_ebeta = self.d_ebeta/1.0e3

    def convert_amu(self):
        """
        Convert micro-amu to amu
        """
        if self.mass:
            self.mass = self.mass/1.0e6
        if self.d_mass:
            self.d_mass = self.d_mass/1.0e6

class AME2012(object):
    def __init__(self, datfile=None):
        """
        Initialize
        """
        self.datfile = datfile
        self.nuclides = []
        
        if self.datfile:
            self.read()

    def read(self):
        """
        Read the AME data file 'datfile'
        """
        try:
            f = open(self.datfile,'r')
        except:
            print('ERROR: data file not found!')
            exit()

        # Get rid of the 39-line header
        for x in range(39):
            f.readline()

        # Read nuclide mass data
        for line in f:
            ls = line.strip()
            if ls:
                # Pull data out of the line
                line, junk = str_head_pop(line, 4)
                line, n = str_head_pop(line, 5)
                line, z = str_head_pop(line, 5)
                line, a = str_head_pop(line, 5)
                line, junk = str_head_pop(line, 1)
                line, element = str_head_pop(line, 3)
                line, origin = str_head_pop(line, 4)
                line, junk = str_head_pop(line, 1)
                line, mexcess = str_head_pop(line, 13)
                mexcess = cleanup(mexcess)
                line, d_mexcess = str_head_pop(line, 11)
                d_mexcess = cleanup(d_mexcess)
                line, nucbind = str_head_pop(line, 11)
                nucbind = cleanup(nucbind)
                line, d_nucbind = str_head_pop(line, 9)
                d_nucbind = cleanup(d_nucbind)
                line, junk = str_head_pop(line, 1)
                line, decay_type = str_head_pop(line, 2)
                line, ebeta = str_head_pop(line, 11)
                ebeta = cleanup(ebeta)
                line, d_ebeta = str_head_pop(line, 9)
                d_ebeta = cleanup(d_ebeta)
                line, junk = str_head_pop(line, 1)
                line, mass = str_head_pop(line, 16)
                mass = cleanup(mass)
                line, d_mass = str_head_pop(line, 11)
                d_mass = cleanup(d_mass)
                # Store nuclide data
                nuclide = Nuclide(n, z, a, element.strip(), origin.strip(),
                                  mexcess, d_mexcess, nucbind,
                                  d_nucbind, decay_type.strip(), ebeta,
                                  d_ebeta, mass, d_mass)
                self.nuclides.append(nuclide)
        f.close()

    def get_nuclide(self, n=-1, z=-1, a=-1):
        """
        Returns the nuclide object given at least 2 of n, z, a.
        
        Complains if A-Z-N != 0.
        """
        if n>=0 and z>=0 and a>=0:
            if a - z - n != 0:
                print('ERROR: Z+N != A')
                exit()
        else:
            if n==-1:
                n = a - z
            elif z==-1:
                z = a - n
            else:
                a = n + z

        for nuc in self.nuclides:
            if (nuc.n == n) and (nuc.z == z):
                return nuc

        print('Nuclide not found for (N, Z, A) = ({}, {}, {})!'.format(n, z, a))
        exit()
