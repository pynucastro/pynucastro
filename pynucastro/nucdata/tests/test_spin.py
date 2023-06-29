from pynucastro.nucdata import SpinTable


class TestSpin:

    @classmethod
    def setup_class(cls):
        """ this is run once for each class before any tests """

    @classmethod
    def teardown_class(cls):
        """ this is run once for each class before any tests """

    def setup_method(self):
        """ this is run once for each class before any tests """
        #Is import to remark that the spin state gs is unique.
        #However, due to some experimental uncertainties two values
        #are proposed.

        self.spintable_gs_reliable = SpinTable(reliable=True)
        self.spintable_gs_not_reliable = SpinTable(reliable=False)

    def teardown_method(self):
        """ this is run once for each class before any tests """

    def test_spin_table(self):

        assert self.spintable_gs_reliable.get_spin_states(a=1, z=0) == 2
        assert self.spintable_gs_reliable.get_spin_states(a=14, z=4) == 1
        assert self.spintable_gs_reliable.get_spin_states(a=219, z=89) == 10
        assert self.spintable_gs_not_reliable.get_spin_states(a=1, z=0) == 2
        assert self.spintable_gs_not_reliable.get_spin_states(a=14, z=4) == 1
        assert self.spintable_gs_not_reliable.get_spin_states(a=219, z=89) == 10
        assert self.spintable_gs_not_reliable.get_spin_states(a=91, z=46) == 8
        assert self.spintable_gs_not_reliable.get_spin_states(a=275, z=107) == 6
        assert self.spintable_gs_not_reliable.get_spin_states(a=11, z=8) == 4
