"""A Fortran wrapper around the SimpleCxxNetwork"""


from pynucastro.networks.simple_cxx_network import SimpleCxxNetwork


class FortranNetwork(SimpleCxxNetwork):
    def __init__(self, *args, **kwargs):

        # Initialize SimpleCxxNetwork parent class
        super().__init__(*args, **kwargs)

    def _get_template_files(self):

        path = self.pynucastro_dir/"templates/fortran-network"

        return path.glob("*.template")
