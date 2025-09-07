import pytest
from pytest import approx

import pynucastro as pyna

# these are some tests that compare directly to the tables
# in the Oda paper.


class TestPruetFullerLibibrary:
    @pytest.fixture(scope="class")
    def pl(self):
        return pyna.PruetFullerLibrary()

    def test_fe66_mn66(self, pl):

        # get the rate
        r = pl.get_rate_by_name("fe66(,)mn66")

        assert r.reactants == [pyna.Nucleus("fe66")]
        assert r.products == [pyna.Nucleus("mn66")]

        # compare to data from the original Pruet & Fuller table,
        # interpolate on a tabulated value, so the interpolation
        # should just be a no-op

        # we need a composition just to define a Ye
        comp = pyna.Composition(["he4", "mn66", "fe66"])
        comp.set_equal()
        Ye = comp.ye

        # from the table -- everything should be -100
        T = 5.e8
        rhoYe = 1.e5
        lbetap = -100
        lec = -100
        lnu = -100

        source_rate = 10.0**-100.0  # 10.0**lbetap + 10.0**lec
        source_nu_loss = 10.0**lnu  # * pyna.constants.constants.MeV2erg

        rate = r.eval(T, rho=rhoYe/Ye, comp=comp)
        assert rate == approx(source_rate, rel=1.e-5, abs=1.e-200)

        nu_loss = r.get_nu_loss(T, rho=rhoYe/Ye, comp=comp)
        assert nu_loss == approx(source_nu_loss, rel=1.e-5, abs=1.e-200)

        # now higher rho and T -- this is not floored to -100
        T = 5.e9
        rhoYe = 1.e8
        lbetap = -16.887
        lec = -13.504
        lnu = -13.310

        source_rate = 10.0**lbetap + 10.0**lec
        source_nu_loss = 10.0**lnu * pyna.constants.constants.MeV2erg

        rate = r.eval(T, rho=rhoYe/Ye, comp=comp)
        assert rate == approx(source_rate, rel=1.e-4, abs=1.e-50)

        nu_loss = r.get_nu_loss(T, rho=rhoYe/Ye, comp=comp)
        assert nu_loss == approx(source_nu_loss, rel=1.e-5, abs=1.e-50)

    def test_fe66_co66(self, pl):

        # get the rate
        r = pl.get_rate_by_name("fe66(,)co66")

        assert r.reactants == [pyna.Nucleus("fe66")]
        assert r.products == [pyna.Nucleus("co66")]

        # compare to data from the original Pruet & Fuller table,
        # interpolate on a tabulated value, so the interpolation
        # should just be a no-op

        # we need a composition just to define a Ye
        comp = pyna.Composition(["he4", "fe66", "co66"])
        comp.set_equal()
        Ye = comp.ye

        # from the table -- everything should be -100
        T = 3.e10
        rhoYe = 1.e9
        lbetam = 0.729
        lpc = 1.466  # positron capture
        lnu = 2.619

        source_rate = 10.0**lbetam + 10.0**lpc
        source_nu_loss = 10.0**lnu * pyna.constants.constants.MeV2erg

        rate = r.eval(T, rho=rhoYe/Ye, comp=comp)
        assert rate == approx(source_rate, rel=1.e-5, abs=1.e-50)

        nu_loss = r.get_nu_loss(T, rho=rhoYe/Ye, comp=comp)
        assert nu_loss == approx(source_nu_loss, rel=1.e-5, abs=1.e-50)

    def test_cu68_zn68(self, pl):

        # get the rate
        r = pl.get_rate_by_name("cu68(,)zn68")

        assert r.reactants == [pyna.Nucleus("cu68")]
        assert r.products == [pyna.Nucleus("zn68")]

        # compare to data from the original Pruet & Fuller table,
        # interpolate on a tabulated value, so the interpolation
        # should just be a no-op

        # we need a composition just to define a Ye
        comp = pyna.Composition(["he4", "cu68", "zn68"])
        comp.set_equal()
        Ye = comp.ye

        # from the table -- everything should be -100
        T = 3.e10
        rhoYe = 1.e11
        lbetam = -2.782
        lpc = -2.423  # positron capture
        lnu = -1.261

        source_rate = 10.0**lbetam + 10.0**lpc
        source_nu_loss = 10.0**lnu * pyna.constants.constants.MeV2erg

        rate = r.eval(T, rho=rhoYe/Ye, comp=comp)
        assert rate == approx(source_rate, rel=1.e-5, abs=1.e-50)

        nu_loss = r.get_nu_loss(T, rho=rhoYe/Ye, comp=comp)
        assert nu_loss == approx(source_nu_loss, rel=1.e-5, abs=1.e-50)
