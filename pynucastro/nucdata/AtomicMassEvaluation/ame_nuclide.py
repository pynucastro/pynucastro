class AMENuclide(object):
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

    def __str__(self):
        return "{}-{}".format(self.element, self.a)
