import pytest
from pytest import approx

import pynucastro as pyna


class TestFFNLibibrary:
    @pytest.fixture(scope="class")
    def ffn_lib(self):
        return pyna.FFNLibrary()

    def test_mg21_na21(self, ffn_lib):

        # get the rate
        r = ffn_lib.get_rate_by_name("mg21(,)na21")

        assert r.reactants == [pyna.Nucleus("mg21")]
        assert r.products == [pyna.Nucleus("na21")]

        # compare to data from the original source (ffnrates.dat) we
        # will interpolate on a tabulated value, so the interpolation
        # should just be a no-op

        # from the table
        T = 0.7e9
        rhoYe = 1.e8
        lbetap = 0.744  # units: log(1/s)
        lepsm = 0.692   # units: log(1/s)
        lnu = 1.829     # units: log(MeV/s)

        source_rate = 10.0**lbetap + 10.0**lepsm
        source_nu_loss = 10.0**lnu * pyna.constants.constants.MeV2erg

        # we need a composition just to define a Ye
        comp = pyna.Composition(["he4", "na21", "mg21"])
        comp.set_equal()
        Ye = comp.ye

        # now interpolate
        rate = r.eval(T, rho=rhoYe/Ye, comp=comp)
        assert rate == approx(source_rate, rel=1.e-4, abs=1.e-30)

        nu_loss = r.get_nu_loss(T, rho=rhoYe/Ye, comp=comp)
        assert nu_loss == approx(source_nu_loss)

    def test_na21_mg21(self, ffn_lib):

        # get the rate
        r = ffn_lib.get_rate_by_name("na21(,)mg21")

        assert r.reactants == [pyna.Nucleus("na21")]
        assert r.products == [pyna.Nucleus("mg21")]

        # compare to data from the original source (ffnrates.dat) we
        # will interpolate on a tabulated value, so the interpolation
        # should just be a no-op

        # from the table
        T = 3e9
        rhoYe = 1.e5
        lbetam = -33.550  # units: log(1/s)
        lepsp = -22.972   # units: log(1/s)
        lnubar = -23.044  # units: log(MeV/s)

        source_rate = 10.0**lbetam + 10.0**lepsp
        source_nu_loss = 10.0**lnubar * pyna.constants.constants.MeV2erg

        # we need a composition just to define a Ye
        comp = pyna.Composition(["he4", "na21", "mg21"])
        comp.set_equal()
        Ye = comp.ye

        # now interpolate
        rate = r.eval(T, rho=rhoYe/Ye, comp=comp)
        assert rate == approx(source_rate, rel=1.e-4, abs=1.e-30)

        nu_loss = r.get_nu_loss(T, rho=rhoYe/Ye, comp=comp)
        assert nu_loss == approx(source_nu_loss)

    def test_p27_si27(self, ffn_lib):

        # get the rate
        r = ffn_lib.get_rate_by_name("p27(,)si27")

        assert r.reactants == [pyna.Nucleus("p27")]
        assert r.products == [pyna.Nucleus("si27")]

        # compare to data from the original source (ffnrates.dat) we
        # will interpolate on a tabulated value, so the interpolation
        # should just be a no-op

        # from the table
        T = 3.e10
        rhoYe = 10.0
        lbetap = 2.259  # units: log(1/s)
        lepsm = 2.052   # units: log(1/s)
        lnu = 3.671     # units: log(MeV/s)

        source_rate = 10.0**lbetap + 10.0**lepsm
        source_nu_loss = 10.0**lnu * pyna.constants.constants.MeV2erg

        # we need a composition just to define a Ye
        comp = pyna.Composition(["he4", "p27", "si27"])
        comp.set_equal()
        Ye = comp.ye

        # now interpolate
        rate = r.eval(T, rho=rhoYe/Ye, comp=comp)
        assert rate == approx(source_rate, rel=1.e-4, abs=1.e-30)

        nu_loss = r.get_nu_loss(T, rho=rhoYe/Ye, comp=comp)
        assert nu_loss == approx(source_nu_loss, rel=1.e-4, abs=1.e-30)

    def test_si27_p27(self, ffn_lib):

        # get the rate
        r = ffn_lib.get_rate_by_name("si27(,)p27")

        assert r.reactants == [pyna.Nucleus("si27")]
        assert r.products == [pyna.Nucleus("p27")]

        # compare to data from the original source (ffnrates.dat) we
        # will interpolate on a tabulated value, so the interpolation
        # should just be a no-op

        # from the table
        T = 3.e10
        rhoYe = 10.0
        lbetam = -1.353  # units: log(1/s)
        lepsp = -0.161   # units: log(1/s)
        lnubar = 0.906     # units: log(MeV/s)

        source_rate = 10.0**lbetam + 10.0**lepsp
        source_nu_loss = 10.0**lnubar * pyna.constants.constants.MeV2erg

        # we need a composition just to define a Ye
        comp = pyna.Composition(["he4", "p27", "si27"])
        comp.set_equal()
        Ye = comp.ye

        # now interpolate
        rate = r.eval(T, rho=rhoYe/Ye, comp=comp)
        assert rate == approx(source_rate, rel=1.e-4, abs=1.e-12)

        nu_loss = r.get_nu_loss(T, rho=rhoYe/Ye, comp=comp)
        assert nu_loss == approx(source_nu_loss, rel=1.e-5, abs=1.e-12)
