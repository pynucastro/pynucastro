from pynucastro.nucdata import MassNuclide, MassTable
from pytest import approx


class TestMass:

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
        self.nuc1 = MassNuclide(a=1, z=0, dm=8.07132)
        self.nuc2 = MassNuclide(a=295, z=118, dm=201.37)

    def teardown_method(self):
        """ this is run once for each class before any tests """
        pass

    def test_nuclide(self):

        assert self.nuc1.a == 1
        assert self.nuc1.z == 0
        assert self.nuc1.dm == approx(8.07132)
        assert str(self.nuc1) == 'n'

        assert self.nuc2.a == 295
        assert self.nuc2.z == 118
        assert self.nuc2.dm == approx(201.37)
        assert str(self.nuc2) == 'og295'

    def test_mass_table(self):

        self._dm = MassTable()
        self.n = self._dm.get_mass_diff('n')
        self.he4 = self._dm.get_mass_diff('he4')
        self.og295 = self._dm.get_mass_diff('og295')

        assert self.nuc1 == self.n
        assert self.nuc2 == self.og295

        assert self.n.a == 1
        assert self.n.z == 0
        assert self.n.dm == approx(8.07132)

        assert self.he4.a == 4
        assert self.he4.z == 2
        assert self.he4.dm == approx(2.42492)

        assert self.og295.a == 295
        assert self.og295.z == 118
        assert self.og295.dm == approx(201.37)
