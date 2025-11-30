import pytest
from pytest import approx

import pynucastro as pyna
from pynucastro import networks
from pynucastro.nucdata import Nucleus

class TestLoddersComposition:

    @pytest.fixture(scope="class")
    def solar(self):
        """Default Lodders solar comp, without scaling"""
        return networks.Lodders()

    @pytest.fixture(scope="class")
    def scaled(self):
        """Scale Lodders to desired metallicity"""
        return networks.Lodders(z=0.04)

    def test_sum_to_one(self, solar, scaled):
        """All mass fractions should add up to 1"""
        xsum_solar = sum(solar.values())
        xsum_scaled = sum(scaled.values())

        assert xsum_solar == approx(1.0)
        assert xsum_scaled == approx(1.0)

    def test_default_z(self, solar):
        """Default metallicity is 1-H-He from Lodders"""

        H = 0.0
        He = 0.0
        metals = 0.0

        for nuc, X in solar.items():
            if nuc.Z==1:
                H += X
            elif nuc.Z==2:
                He += X
            elif nuc.Z>=3:
                metals += X

        assert metals == approx(0.014964158698946859)

    def test_scaled_comp(self, scaled):
        """Checking for the test value Z=0.04"""

        scaled_H = 0.0
        scaled_He = 0.0
        scaled_z = 0.0

        for nuc, X in scaled.items():
            if nuc.Z == 1:
                scaled_H += X
            elif nuc.Z == 2:
                scaled_He += X
            elif nuc.Z >= 3:
                scaled_z += X

        assert scaled_H == approx( 0.7271232554513777)
        assert scaled_He == approx(0.23276984006179305)
        assert scaled_z == approx(0.0401069044868292)





