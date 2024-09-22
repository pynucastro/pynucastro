from pytest import approx

from pynucastro.nucdata import MassTable


class TestMass:
    def test_mass_table(self):
        _dm = MassTable()
        assert _dm.get_mass_diff(a=1, z=0) == approx(8.0713181)
        assert _dm.get_mass_diff(a=4, z=2) == approx(2.42491587)
        assert _dm.get_mass_diff(a=295, z=118) == approx(201.37)
