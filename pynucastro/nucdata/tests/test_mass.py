from pynucastro.nucdata import MassNuclide, MassTable
import pytest
from pytest import approx


class TestMass:
    @pytest.fixture(scope="class")
    def nuc1(self):
        return MassNuclide(a=1, z=0, dm=8.07132)

    @pytest.fixture(scope="class")
    def nuc2(self):
        return MassNuclide(a=295, z=118, dm=201.37)

    def test_nuclide(self, nuc1, nuc2):

        assert nuc1.a == 1
        assert nuc1.z == 0
        assert nuc1.dm == approx(8.07132)
        assert str(nuc1) == 'n'

        assert nuc2.a == 295
        assert nuc2.z == 118
        assert nuc2.dm == approx(201.37)
        assert str(nuc2) == 'og295'

    def test_mass_table(self, nuc1, nuc2):

        _dm = MassTable()
        n = _dm.get_mass_diff('n')
        he4 = _dm.get_mass_diff('he4')
        og295 = _dm.get_mass_diff('og295')

        assert nuc1 == n
        assert nuc2 == og295

        assert n.a == 1
        assert n.z == 0
        assert n.dm == approx(8.07132)

        assert he4.a == 4
        assert he4.z == 2
        assert he4.dm == approx(2.42492)

        assert og295.a == 295
        assert og295.z == 118
        assert og295.dm == approx(201.37)
