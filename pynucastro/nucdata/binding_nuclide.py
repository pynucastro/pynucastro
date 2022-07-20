# Common Imports

class BindingNuclide:
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
        print(f'n = {self.n}')
        print(f'z = {self.z}')
        print(f'nucbind = {self.nucbind}')
