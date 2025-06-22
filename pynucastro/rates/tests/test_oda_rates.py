import pytest
from pytest import approx

import pynucastro as pyna

# these are some tests that compare directly to the tables
# in the Oda paper.


class TestOdaLibibrary:
    @pytest.fixture(scope="class")
    def oda_lib(self):
        return pyna.OdaLibrary()

    def test_f17_o17(self, oda_lib):

        # get the rate
        r = oda_lib.get_rate_by_name("f17(,)o17")

        assert r.reactants == [pyna.Nucleus("f17")]
        assert r.products == [pyna.Nucleus("o17")]

        # compare to data from the original Oda table will interpolate
        # on a tabulated value, so the interpolation should just be a
        # no-op

        # from the table
        T = 0.7e9
        rhoYe = 1.e8
        lbetap = -1.986  # units: log(1/s)
        lepsm = -0.079   # units: log(1/s)
        lnu = 0.547     # units: log(MeV/s)

        source_rate = 10.0**lbetap + 10.0**lepsm
        source_nu_loss = 10.0**lnu * pyna.constants.constants.MeV2erg

        # we need a composition just to define a Ye
        comp = pyna.Composition(["he4", "o17", "f17"])
        comp.set_equal()
        Ye = comp.ye

        # now interpolate
        rate = r.eval(T, rho=rhoYe/Ye, comp=comp)
        assert rate == approx(source_rate, rel=1.e-4, abs=1.e-30)

        nu_loss = r.get_nu_loss(T, rho=rhoYe/Ye, comp=comp)
        assert nu_loss == approx(source_nu_loss)

    def test_o17_f17(self, oda_lib):

        # get the rate
        r = oda_lib.get_rate_by_name("o17(,)f17")

        assert r.reactants == [pyna.Nucleus("o17")]
        assert r.products == [pyna.Nucleus("f17")]

        # compare to data from the original Oda table will interpolate
        # on a tabulated value, so the interpolation should just be a
        # no-op

        # from the table
        T = 3.e9
        rhoYe = 1.e9
        lbetam = -15.286  # units: log(1/s)
        lepsp = -15.041   # units: log(1/s)
        lnu = -14.933     # units: log(MeV/s)

        source_rate = 10.0**lbetam + 10.0**lepsp
        source_nu_loss = 10.0**lnu * pyna.constants.constants.MeV2erg

        # we need a composition just to define a Ye
        comp = pyna.Composition(["he4", "o17", "f17"])
        comp.set_equal()
        Ye = comp.ye

        # now interpolate
        rate = r.eval(T, rho=rhoYe/Ye, comp=comp)
        assert rate == approx(source_rate, rel=1.e-4, abs=1.e-30)

        nu_loss = r.get_nu_loss(T, rho=rhoYe/Ye, comp=comp)
        assert nu_loss == approx(source_nu_loss)

    def test_al23_mg23(self, oda_lib):

        # get the rate
        r = oda_lib.get_rate_by_name("al23(,)mg23")

        assert r.reactants == [pyna.Nucleus("al23")]
        assert r.products == [pyna.Nucleus("mg23")]

        # compare to data from the original Oda table will interpolate
        # on a tabulated value, so the interpolation should just be a
        # no-op

        # from the table
        T = 1.e10
        rhoYe = 1.e8
        lbetap = 0.546  # units: log(1/s)
        lepsm = 0.902   # units: log(1/s)
        lnu = 1.929     # units: log(MeV/s)

        source_rate = 10.0**lbetap + 10.0**lepsm
        source_nu_loss = 10.0**lnu * pyna.constants.constants.MeV2erg

        # we need a composition just to define a Ye
        comp = pyna.Composition(["he4", "al23", "mg23"])
        comp.set_equal()
        Ye = comp.ye

        # now interpolate
        rate = r.eval(T, rho=rhoYe/Ye, comp=comp)
        assert rate == approx(source_rate, rel=1.e-4, abs=1.e-30)

        nu_loss = r.get_nu_loss(T, rho=rhoYe/Ye, comp=comp)
        assert nu_loss == approx(source_nu_loss)

    def test_mg23_al23(self, oda_lib):

        # get the rate
        r = oda_lib.get_rate_by_name("mg23(,)al23")

        assert r.reactants == [pyna.Nucleus("mg23")]
        assert r.products == [pyna.Nucleus("al23")]

        # compare to data from the original Oda table will interpolate
        # on a tabulated value, so the interpolation should just be a
        # no-op

        # from the table
        T = 1.e10
        rhoYe = 1.e8
        lbetam = -8.609  # units: log(1/s)
        lepsp = -6.611   # units: log(1/s)
        lnu = -6.082     # units: log(MeV/s)

        source_rate = 10.0**lbetam + 10.0**lepsp
        source_nu_loss = 10.0**lnu * pyna.constants.constants.MeV2erg

        # we need a composition just to define a Ye
        comp = pyna.Composition(["he4", "mg23", "al23"])
        comp.set_equal()
        Ye = comp.ye

        # now interpolate
        rate = r.eval(T, rho=rhoYe/Ye, comp=comp)
        assert rate == approx(source_rate, rel=1.e-4, abs=1.e-30)

        nu_loss = r.get_nu_loss(T, rho=rhoYe/Ye, comp=comp)
        assert nu_loss == approx(source_nu_loss)

    def test_p31_si31(self, oda_lib):

        # get the rate
        r = oda_lib.get_rate_by_name("p31(,)si31")

        assert r.reactants == [pyna.Nucleus("p31")]
        assert r.products == [pyna.Nucleus("si31")]

        # compare to data from the original Oda table will interpolate
        # on a tabulated value, so the interpolation should just be a
        # no-op

        # from the table
        T = 5.e9
        rhoYe = 1.e5
        lbetap = -6.515  # units: log(1/s)
        lepsm = -4.850   # units: log(1/s)
        lnu = -4.528     # units: log(MeV/s)

        source_rate = 10.0**lbetap + 10.0**lepsm
        source_nu_loss = 10.0**lnu * pyna.constants.constants.MeV2erg

        # we need a composition just to define a Ye
        comp = pyna.Composition(["he4", "p31", "si31"])
        comp.set_equal()
        Ye = comp.ye

        # now interpolate
        rate = r.eval(T, rho=rhoYe/Ye, comp=comp)
        assert rate == approx(source_rate, rel=1.e-4, abs=1.e-30)

        nu_loss = r.get_nu_loss(T, rho=rhoYe/Ye, comp=comp)
        assert nu_loss == approx(source_nu_loss)

    def test_si31_p31(self, oda_lib):

        # get the rate
        r = oda_lib.get_rate_by_name("si31(,)p31")

        assert r.reactants == [pyna.Nucleus("si31")]
        assert r.products == [pyna.Nucleus("p31")]

        # compare to data from the original Oda table will interpolate
        # on a tabulated value, so the interpolation should just be a
        # no-op

        # from the table
        T = 5.e9
        rhoYe = 1.e5
        lbetam = -3.412  # units: log(1/s)
        lepsp = -3.478   # units: log(1/s)
        lnu = -2.817     # units: log(MeV/s)

        source_rate = 10.0**lbetam + 10.0**lepsp
        source_nu_loss = 10.0**lnu * pyna.constants.constants.MeV2erg

        # we need a composition just to define a Ye
        comp = pyna.Composition(["he4", "si31", "p31"])
        comp.set_equal()
        Ye = comp.ye

        # now interpolate
        rate = r.eval(T, rho=rhoYe/Ye, comp=comp)
        assert rate == approx(source_rate, rel=1.e-4, abs=1.e-30)

        nu_loss = r.get_nu_loss(T, rho=rhoYe/Ye, comp=comp)
        assert nu_loss == approx(source_nu_loss)
