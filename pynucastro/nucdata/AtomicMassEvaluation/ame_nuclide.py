# Common Imports

class AMENuclide:
    """Hold the information for a single nucleus from the AME
    database."""
    def __init__(self, n=None, z=None, a=None, element=None, origin=None,
                 mexcess=None, d_mexcess=None, nucbind=None,
                 d_nucbind=None, decay_type=None, ebeta=None,
                 d_ebeta=None, mass=None, d_mass=None):

        self.n = None
        self.z = None
        self.a = None
        self.element = None
        self.origin = None
        self.mexcess = None
        self.d_mexcess = None
        self.nucbind = None
        self.d_nucbind = None
        self.decay_type = None
        self.ebeta = None
        self.d_ebeta = None
        self.mass = None
        self.d_mass = None

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
        print(f'n = {self.n}')
        print(f'z = {self.z}')
        print(f'a = {self.a}')
        print(f'element = {self.element}')
        print(f'origin = {self.origin}')
        print(f'mexcess = {self.mexcess}')
        print(f'd_mexcess = {self.d_mexcess}')
        print(f'nucbind = {self.nucbind}')
        print(f'd_nucbind = {self.d_nucbind}')
        print(f'decay_type = {self.decay_type}')
        print(f'ebeta = {self.ebeta}')
        print(f'd_ebeta = {self.d_ebeta}')
        print(f'mass = {self.mass}')
        print(f'd_mass = {self.d_mass}')

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

    def __str__(self):
        return f"{self.element}-{self.a}"
