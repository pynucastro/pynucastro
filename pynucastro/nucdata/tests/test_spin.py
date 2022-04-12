from pynucastro.nucdata import SpinNuclide, SpinTable
from numpy import array


class TestSpin:

    @classmethod
    def setup_class(cls):
        """ this is run once for each class before any tests """
        pass

    @classmethod
    def teardown_class(cls):
        """ this is run once for each class before any tests """
        pass

    def setup_method(self):
        """ this is run once for each class before any tests """
        #Is import to remark that the spin state gs is unique.
        #However, due to some experimental uncertainties two values
        #are proposed.

        self.spintable_gs = SpinTable(set_double_gs=False)
        self.spintable_double_gs = SpinTable(set_double_gs=True)

    def teardown_method(self):
        """ this is run once for each class before any tests """
        pass

    def test_spin_table(self):

        #assert self.spintable_gs.get_spin_nuclide('n') == SpinNuclide(a=1, z=0, spin_states=2)
        assert self.spintable_gs.get_spin_nuclide('pd91') == SpinNuclide(a=91, z=46, spin_states=8)
        assert self.spintable_double_gs.get_spin_nuclide('br89') == SpinNuclide(a=89, z=35, spin_states=array([4, 6]))
        assert self.spintable_double_gs.get_spin_nuclide('pr130') == SpinNuclide(a=130, z=59, spin_states=array([13, 15]))
