"""A Fortran wrapper around the SimpleCxxNetwork"""


from pynucastro.networks.simple_cxx_network import SimpleCxxNetwork


class FortranNetwork(SimpleCxxNetwork):
    def __init__(self, *args, **kwargs):

        # Initialize SimpleCxxNetwork parent class
        super().__init__(*args, **kwargs)

        self.ftags['<nspec_fortran>'] = self._nspec_fortran
        self.ftags['<aion_fortran>'] = self._aion_fortran
        self.ftags['<zion_fortran>'] = self._zion_fortran
        self.ftags['<name_fortran>'] = self._name_fortran

    def _get_template_files(self):

        path = self.pynucastro_dir/"templates/fortran-network"

        return path.glob("*.template")

    def _nspec_fortran(self, n_indent, of):
        of.write(f'{self.indent*n_indent}integer, parameter :: nspec = {len(self.unique_nuclei)}\n')

    def _aion_fortran(self, n_indent, of):
        of.write(f'{self.indent*n_indent}double precision, parameter :: aion(nspec) = [ &\n')
        for n, nuc in enumerate(self.unique_nuclei):
            if n % 4 == 0:
                of.write('        ')
            of.write(f'{nuc.A}')
            if n < len(self.unique_nuclei)-1:
                of.write(', ')
            if n % 4 == 3:
                of.write(' &\n')
        if len(self.unique_nuclei) % 4 == 0:
            of.write('       ')
        of.write(']\n')

    def _zion_fortran(self, n_indent, of):
        of.write(f'{self.indent*n_indent}double precision, parameter :: zion(nspec) = [ &\n')
        for n, nuc in enumerate(self.unique_nuclei):
            if n % 4 == 0:
                of.write('        ')
            of.write(f'{nuc.Z}')
            if n < len(self.unique_nuclei)-1:
                of.write(', ')
            if n % 4 == 3:
                of.write(' &\n')
        if len(self.unique_nuclei) % 4 == 0:
            of.write('       ')
        of.write(']\n')

    def _name_fortran(self, n_indent, of):
        of.write(f'{self.indent*n_indent}character(len=6), dimension(nspec) :: spec_names = [ &\n')
        for n, nuc in enumerate(self.unique_nuclei):
            if n == 0:
                of.write('        character(len=6) :: ')
            elif n % 4 == 0:
                of.write('                            ')
            of.write(f'"{nuc.short_spec_name}"')
            if n < len(self.unique_nuclei)-1:
                of.write(', ')
            if n % 4 == 3:
                of.write(' &\n')
        if len(self.unique_nuclei) % 4 == 0:
            of.write('       ')
        of.write(']\n')
