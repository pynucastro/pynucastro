from pynucastro.nucdata import SpinNuclide, SpinTable


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

        self.spintable_gs_reliable = SpinTable(reliable=True)
        self.spintable_gs_not_reliable = SpinTable(reliable=False)

    def teardown_method(self):
        """ this is run once for each class before any tests """
        pass

    def test_spin_table(self):

        assert self.spintable_gs_reliable.get_spin_nuclide('n') == SpinNuclide(a=1, z=0, spin_states=2)
        assert self.spintable_gs_reliable.get_spin_nuclide('be14') == SpinNuclide(a=14, z=4, spin_states=1)
        assert self.spintable_gs_reliable.get_spin_nuclide('ac219') == SpinNuclide(a=219, z=89, spin_states=10)
        assert self.spintable_gs_not_reliable.get_spin_nuclide('pd91') == SpinNuclide(a=91, z=46, spin_states=8)
        assert self.spintable_gs_not_reliable.get_spin_nuclide('bh275') == SpinNuclide(a=275, z=107, spin_states=6)
        assert self.spintable_gs_not_reliable.get_spin_nuclide('o11') == SpinNuclide(a=11, z=8, spin_states=4)
