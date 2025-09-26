import pytest

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

        self.spintable_gs = SpinTable()

    def teardown_method(self):
        """ this is run once for each class before any tests """

    def test_spin_table(self):
        #Is important to remark that the spin state gs is unique.
        #However, due to some experimental uncertainties two values
        #are proposed.

        assert self.spintable_gs.get_spin_states(a=1, z=0, reliable=True) == 2
        assert self.spintable_gs.get_spin_states(a=14, z=4, reliable=True) == 1
        assert self.spintable_gs.get_spin_states(a=219, z=89, reliable=True) == 10
        with pytest.raises(NotImplementedError):
            self.spintable_gs.get_spin_states(a=91, z=46, reliable=True)
        with pytest.raises(NotImplementedError):
            self.spintable_gs.get_spin_states(a=275, z=107, reliable=True)
        with pytest.raises(NotImplementedError):
            self.spintable_gs.get_spin_states(a=11, z=8, reliable=True)

        assert self.spintable_gs.get_spin_states(a=1, z=0, reliable=False) == 2
        assert self.spintable_gs.get_spin_states(a=14, z=4, reliable=False) == 1
        assert self.spintable_gs.get_spin_states(a=219, z=89, reliable=False) == 10
        assert self.spintable_gs.get_spin_states(a=91, z=46, reliable=False) == 8
        assert self.spintable_gs.get_spin_states(a=275, z=107, reliable=False) == 6
        assert self.spintable_gs.get_spin_states(a=11, z=8, reliable=False) == 4
