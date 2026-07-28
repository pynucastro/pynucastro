# unit tests for different types of rates

import copy
import math

import pytest
from pytest import approx

from pynucastro import Composition, Rate, ReacLibLibrary, rates
from pynucastro.nucdata import Nucleus
from pynucastro.rates import BaryonConservationError, ModifiedRate
from pynucastro.rates.alternate_rates import IliadisO16pgF17
from pynucastro.screening.screen import chugunov_2007


class TestTfactors:
    @pytest.fixture(scope="class")
    @classmethod
    def tf(cls):
        return rates.Tfactors(2.e9)

    def test_tfactors(self, tf):
        assert tf.T9 == approx(2.0)
        assert tf.T9i == approx(0.5)
        assert tf.T913i == approx(0.5**(1./3.))
        assert tf.T913 == approx(2.0**(1./3.))
        assert tf.T953 == approx(2.0**(5./3.))
        assert tf.lnT9 == approx(math.log(2.0))


class TestRate:

    @pytest.fixture(scope="class")
    @classmethod
    def srates(cls, reaclib_library):

        sample_rates = {}

        # chapter-1
        sample_rates["ch1"] = reaclib_library.get_rate_by_name("o15(,)n15")

        # chapter-2
        sample_rates["ch2"] = reaclib_library.get_rate_by_name("t(g,n)d")

        # chapter-3
        sample_rates["ch3"] = reaclib_library.get_rate_by_name("he6(g,nn)he4")

        # chapter-4
        sample_rates["ch4"] = reaclib_library.get_rate_by_name("c12(a,g)o16")

        # chapter-5
        sample_rates["ch5"] = reaclib_library.get_rate_by_name("n15(p,a)c12")

        # chapter-6
        sample_rates["ch6"] = reaclib_library.get_rate_by_name("he3(he3,pp)he4")

        # chapter-7
        sample_rates["ch7"] = reaclib_library.get_rate_by_name("li7(t,nna)he4")

        # chapter-8
        sample_rates["ch8"] = reaclib_library.get_rate_by_name("he4(aa,g)c12")

        # chapter-9
        sample_rates["ch9"] = reaclib_library.get_rate_by_name("he4(pp,he3)he3")

        # chapter-10
        sample_rates["ch10"] = reaclib_library.get_rate_by_name("he4(npa,he3)li7")

        # chapter-11
        sample_rates["ch11"] = reaclib_library.get_rate_by_name("b17(g,nnn)c14")

        return sample_rates

    def test_source(self, srates):
        assert srates["ch1"].source["Year"] == "2012"

    def test_chapter(self, srates):
        assert srates["ch1"].chapter == 1
        assert srates["ch2"].chapter == 2
        assert srates["ch3"].chapter == 3
        assert srates["ch4"].chapter == 4
        assert srates["ch5"].chapter == 5
        assert srates["ch6"].chapter == 6
        assert srates["ch7"].chapter == 7
        assert srates["ch8"].chapter == 8
        assert srates["ch9"].chapter == 9
        assert srates["ch10"].chapter == 10
        assert srates["ch11"].chapter == 11

    def test_labelprops(self, srates):
        assert srates["ch1"].labelprops == "wc12w "
        assert srates["ch2"].labelprops == "nk06nv"
        assert srates["ch3"].labelprops == "cf88rv"
        assert srates["ch4"].labelprops == "nac2  "
        assert srates["ch5"].labelprops == "nacrc "
        assert srates["ch6"].labelprops == "nacrn "
        assert srates["ch7"].labelprops == "mafon "
        assert srates["ch8"].labelprops == "fy05c "
        assert srates["ch9"].labelprops == "nacrnv"
        assert srates["ch10"].labelprops == "mafonv"
        assert srates["ch11"].labelprops == "wc12w "

    def test_reactants(self, srates):

        # o15--n15-wc12
        assert srates["ch1"].reactants[0] == Nucleus("o15")
        assert len(srates["ch1"].reactants) == 1

        # t-gn-d-nk06
        assert srates["ch2"].reactants[0] == Nucleus("h3")
        assert len(srates["ch2"].reactants) == 1

        # he6-gnn-he4-cf88
        assert srates["ch3"].reactants[0] == Nucleus("he6")
        assert len(srates["ch3"].reactants) == 1

        # c12-ag-o16-nac2
        assert srates["ch4"].reactants[0] == Nucleus("he4")
        assert srates["ch4"].reactants[1] == Nucleus("c12")
        assert len(srates["ch4"].reactants) == 2

        # n15-pa-c12-nacr
        assert srates["ch5"].reactants[0] == Nucleus("h1")
        assert srates["ch5"].reactants[1] == Nucleus("n15")
        assert len(srates["ch5"].reactants) == 2

        # he3-he3pp-he4-nacr
        assert srates["ch6"].reactants[0] == Nucleus("he3")
        assert srates["ch6"].reactants[1] == Nucleus("he3")
        assert len(srates["ch6"].reactants) == 2

        # li7-tnna-he4-mafo
        assert srates["ch7"].reactants[0] == Nucleus("h3")
        assert srates["ch7"].reactants[1] == Nucleus("li7")
        assert len(srates["ch7"].reactants) == 2

        # he4-aag-c12-fy05
        assert srates["ch8"].reactants[0] == Nucleus("he4")
        assert srates["ch8"].reactants[1] == Nucleus("he4")
        assert srates["ch8"].reactants[2] == Nucleus("he4")
        assert len(srates["ch8"].reactants) == 3

        # he4-pphe3-he3-nacr
        assert srates["ch9"].reactants[0] == Nucleus("p")
        assert srates["ch9"].reactants[1] == Nucleus("h1")
        assert srates["ch9"].reactants[2] == Nucleus("he4")
        assert len(srates["ch9"].reactants) == 3

        # he4-npahe3-li7-mafo
        assert srates["ch10"].reactants[0] == Nucleus("n")
        assert srates["ch10"].reactants[1] == Nucleus("h1")
        assert srates["ch10"].reactants[2] == Nucleus("he4")
        assert srates["ch10"].reactants[3] == Nucleus("he4")
        assert len(srates["ch10"].reactants) == 4

        # b17-nnn-c14-wc12
        assert srates["ch11"].reactants[0] == Nucleus("b17")
        assert len(srates["ch11"].reactants) == 1

    def test_products(self, srates):
        assert srates["ch4"].products[0] == Nucleus("o16")
        assert srates["ch8"].products[0] == Nucleus("c12")
        assert len(srates["ch8"].products) == 1

        # he4-pphe3-he3-nacr
        assert srates["ch9"].products[0] == Nucleus("he3")
        assert srates["ch9"].products[1] == Nucleus("he3")
        assert len(srates["ch9"].products) == 2

    def test_count(self, srates):
        assert srates["ch3"].reactant_count(Nucleus("he6")) == 1
        assert srates["ch3"].product_count(Nucleus("he6")) == 0

        assert srates["ch5"].reactant_count(Nucleus("n15")) == 1
        assert srates["ch5"].reactant_count(Nucleus("p")) == 1
        assert srates["ch5"].reactant_count(Nucleus("a")) == 0
        assert srates["ch5"].reactant_count(Nucleus("c12")) == 0

        assert srates["ch5"].product_count(Nucleus("n15")) == 0
        assert srates["ch5"].product_count(Nucleus("p")) == 0
        assert srates["ch5"].product_count(Nucleus("a")) == 1
        assert srates["ch5"].product_count(Nucleus("c12")) == 1

        assert srates["ch8"].reactant_count(Nucleus("a")) == 3

    def test_prefactor(self, srates):
        assert srates["ch4"].prefactor == 1.0
        assert srates["ch8"].prefactor == approx(0.16666666)

    def test_rate_exponent(self, srates):
        assert srates["ch8"].get_rate_exponent(1.e8) == approx(40.9106396)

    def test_eval(self, srates):
        assert srates["ch8"].eval(1.e8) == approx(2.0403192412842946e-24, rel=1.e-6, abs=1.e-40)

    def test_eval_deriv(self, srates):
        T0 = 1.e8
        eps = 1.e-8

        # compare finite diff to analytic diff

        # rate4
        diff = (srates["ch4"].eval(T0*(1.0+eps)) - srates["ch4"].eval(T0)) / (T0 * eps)
        err = abs(diff - srates["ch4"].eval_deriv(T0)) / diff

        assert err < 1.e-6

        # rate5
        diff = (srates["ch5"].eval(T0*(1.0+eps)) - srates["ch5"].eval(T0)) / (T0 * eps)
        err = abs(diff - srates["ch5"].eval_deriv(T0)) / diff

        assert err < 1.e-6

        # rate6
        diff = (srates["ch6"].eval(T0*(1.0+eps)) - srates["ch6"].eval(T0)) / (T0 * eps)
        err = abs(diff - srates["ch6"].eval_deriv(T0)) / diff

        assert err < 1.e-6

    def test_comparison(self, srates):
        assert srates["ch1"] > srates["ch2"]
        assert srates["ch1"] > srates["ch4"]
        assert srates["ch8"] > srates["ch9"]

    def test_weak(self, srates):
        assert srates["ch1"].weak
        assert not srates["ch2"].weak

    def test_weak_type(self, srates):
        assert srates["ch1"].weak_type == "beta_pos"
        assert srates["ch2"].weak_type == ""
        assert srates["ch11"].weak_type == "beta_neg"

    def test_screen(self, srates):
        assert not srates["ch1"].ion_screen
        assert srates["ch4"].ion_screen == [Nucleus("he4"), Nucleus("c12")]
        assert srates["ch8"].ion_screen == 3*[Nucleus("he4")]

    def test_heaviest_lightest(self, srates):
        assert srates["ch4"].heaviest() == Nucleus("o16")
        assert srates["ch4"].lightest() == Nucleus("he4")
        assert srates["ch2"].lightest() == Nucleus("n")
        assert srates["ch2"].heaviest() == Nucleus("t")

    def test_identical_particle_factor(self, srates):

        # work on a copy so this doesn't persist in the library as a
        # reference

        rr = copy.deepcopy(srates["ch8"])

        assert rr.prefactor == approx(0.16666667)

        rr.use_identical_particle_factor = False
        rr._set_rhs_properties()  # pylint: disable=protected-access

        assert rr.prefactor == 1.0

    def test_stoichiometry(self, srates, reaclib_library):
        assert repr(srates["ch4"]) == "C12 + He4 ⟶ O16 + 𝛾"

        # create a separate version since rates are mutable
        _c12ag = reaclib_library.get_rate_by_name("c12(a,g)o16")
        c12ag = copy.deepcopy(_c12ag)

        c12ag.stoichiometry = {Nucleus("he4"): 1.5,
                               Nucleus("c12"): 1,
                               Nucleus("o16"): 1.125}
        c12ag._set_print_representation()  # pylint: disable=protected-access

        assert repr(c12ag) == "C12 + 1.5 He4 ⟶ 1.125 O16 + 𝛾"
        assert c12ag.rid == "C12 + 1.5 He4 --> 1.125 O16"

        assert c12ag.reactant_count(Nucleus("he4")) == 1.5

    def test_stoichiometry_3alpha(self, reaclib_library):

        _three_alpha = reaclib_library.get_rate_by_name("he4(aa,g)c12")
        three_alpha = copy.deepcopy(_three_alpha)

        assert repr(three_alpha) == "3 He4 ⟶ C12 + 𝛾"
        assert three_alpha.rid == "3 He4 --> C12"

        three_alpha.stoichiometry = {Nucleus("he4"): 4,
                                     Nucleus("c12"): 4/3}
        three_alpha._set_print_representation()  # pylint: disable=protected-access

        assert repr(three_alpha) == "4 He4 ⟶ 1.3333333333333333 C12 + 𝛾"
        assert three_alpha.rid == "4 He4 --> 1.3333333333333333 C12"

        assert three_alpha.reactant_count(Nucleus("he4")) == 4

    def test_stoichiometry_dict(self, reaclib_library):

        _c12c12 = reaclib_library.get_rate_by_name("c12(c12,a)ne20")
        c12c12 = copy.deepcopy(_c12c12)

        c12c12.stoichiometry = {Nucleus("he4"): 4,
                                Nucleus("ne20"): 0.4}

        c12c12._set_print_representation()  # pylint: disable=protected-access

        assert repr(c12c12) == "C12 + C12 ⟶ 4 He4 + 0.4 Ne20"
        assert c12c12.rid == "C12 + C12 --> 4 He4 + 0.4 Ne20"

        assert c12c12.reactant_count(Nucleus("he4")) == 0
        assert c12c12.product_count(Nucleus("he4")) == 4

        assert c12c12.reactant_count(Nucleus("c12")) == 2

        assert c12c12.product_count(Nucleus("ne20")) == 0.4

    def test_baryon_conservation(self):

        # this will raise an exception
        with pytest.raises(BaryonConservationError):
            _ = Rate(reactants=[Nucleus("n14"), Nucleus("he4")],
                     products=[Nucleus("ne20")])

        # this conserves baryon number
        _ = Rate(reactants=[Nucleus("n14"), Nucleus("he4")],
                 products=[Nucleus("ne20")],
                 stoichiometry={Nucleus("he4"): 1.5})


class TestDerivedRate:

    def test_a_a_ag_c12(self, reaclib_library):
        """
        Here we test the inverse rate, computed by the use of detailed balance
        of a:

        A + B -> C

        reaction type.
        """

        a_a_ag_c12 = reaclib_library.get_rate_by_name("he4(aa,g)c12")
        c12_ga_a_a_reaclib = reaclib_library.get_rate_by_name("c12(g,aa)he4")
        c12_ga_a_a_derived = rates.DerivedRate(source_rate=a_a_ag_c12, use_pf=False, use_unreliable_spins=False)

        assert c12_ga_a_a_reaclib.eval(T=2.0e9) == approx(c12_ga_a_a_derived.eval(T=2.0e9), rel=7.e-3)

    def test_a_a_ag_c12_with_pf(self, reaclib_library):
        """
        This function test the correct rate value if we take in consideration the partition
        functions on the range 1.0e9 to 100.0e9
        """

        a_a_ag_c12 = reaclib_library.get_rate_by_name("he4(aa,g)c12")
        c12_ga_a_a_derived = rates.DerivedRate(source_rate=a_a_ag_c12, use_pf=True, use_unreliable_spins=False)

        with pytest.warns(UserWarning, match="C12 partition function is not supported by tables"):
            rval = c12_ga_a_a_derived.eval(T=2.0e9)
        assert rval == approx(2.9138256017033057e-07)

    def test_iliadis_o16_pg_f17_derived(self):
        o16pgf17 = IliadisO16pgF17()
        f17gpo16_derived = rates.DerivedRate(source_rate=o16pgf17, use_pf=True,
                                             use_unreliable_spins=False)

        rval = f17gpo16_derived.eval(T=2.0e9)
        assert rval == approx(20489931419.399235)


class TestWeakRates:
    @pytest.fixture(scope="class")
    @classmethod
    def rate1(cls):
        return rates.TabularRate("suzuki-18o-18f_betadecay.dat")

    @pytest.fixture(scope="class")
    @classmethod
    def rate2(cls):
        return rates.TabularRate("suzuki-22na-22ne_electroncapture.dat")

    @pytest.fixture(scope="class")
    @classmethod
    def rate3(cls):
        return rates.TabularRate("langanke-45sc-45ca_electroncapture.dat")

    @pytest.fixture(scope="class")
    @classmethod
    def rate4(cls):
        return rates.TabularRate("langanke-45ti-45sc_electroncapture.dat")

    @pytest.fixture(scope="class")
    @classmethod
    def rate5(cls):
        return rates.TabularRate("langanke-45v-45ti_electroncapture.dat")

    @pytest.fixture(scope="class")
    @classmethod
    def rate6(cls):
        return rates.TabularRate("langanke-45ca-45sc_betadecay.dat")

    def test_reactants(self, rate1, rate2, rate3, rate4, rate5, rate6):

        # pick a composition that gives Ye = 0.5 just for testing
        comp = Composition(["c12", "o16"])
        comp.set_equal()

        assert len(rate1.reactants) == 1 and len(rate1.products) == 1
        assert rate1.products[0] == Nucleus("f18")
        assert rate1.reactants[0] == Nucleus("o18")
        assert rate1.eval(2.5e9, rho=2.e8, comp=comp) == approx(8.032467196099662e-16, rel=1.e-6, abs=1.e-20)

        assert len(rate2.reactants) == 1 and len(rate2.products) == 1
        assert rate2.products[0] == Nucleus("ne22")
        assert rate2.reactants[0] == Nucleus("na22")
        assert rate2.eval(1.e9, rho=4.e7, comp=comp) == approx(3.232714235735518e-05, rel=1.e-6, abs=1.e-20)

        assert len(rate3.reactants) == 1 and len(rate3.products) == 1
        assert rate3.products[0] == Nucleus("ca45")
        assert rate3.reactants[0] == Nucleus("sc45")
        assert math.log10(rate3.eval(1.e9, rho=2.e11, comp=comp)) == approx(3.4400000000000004)

        assert len(rate4.reactants) == 1 and len(rate4.products) == 1
        assert rate4.products[0] == Nucleus("sc45")
        assert rate4.reactants[0] == Nucleus("ti45")
        assert math.log10(rate4.eval(1.e9, rho=2.e11, comp=comp)) == approx(3.853)

        assert len(rate5.reactants) == 1 and len(rate5.products) == 1
        assert rate5.products[0] == Nucleus("ti45")
        assert rate5.reactants[0] == Nucleus("v45")
        assert math.log10(rate5.eval(1.e9, rho=2.e11, comp=comp)) == approx(4.71501)

        assert len(rate6.reactants) == 1 and len(rate6.products) == 1
        assert rate6.products[0] == Nucleus("sc45")
        assert rate6.reactants[0] == Nucleus("ca45")
        assert math.log10(rate6.eval(1.e9, rho=2.e11, comp=comp)) == approx(-99.69797)


class TestModify:
    @pytest.fixture(scope="function")
    @classmethod
    def rate(cls, reaclib_library):
        return copy.deepcopy(reaclib_library.get_rate_by_name("c12(c12,n)mg23"))

    def test_modify(self, rate):

        rate.modify_products("mg24")

        assert rate.Q == approx(13.933578000000125)
        assert rate.products == [Nucleus("mg24")]
        assert rate.modified


class TestModifiedRate:
    @pytest.fixture(scope="function")
    @classmethod
    def rate(cls, reaclib_library):
        return copy.deepcopy(reaclib_library.get_rate_by_name("c12(c12,n)mg23"))

    def test_eval(self, rate):

        # pick a composition that gives Ye = 0.5 just for testing
        comp = Composition(["c12", "o16"])
        comp.set_equal()

        # Test eval against the original rate without updating ion_screen,
        # i.e. use same reactants as the original rate
        r1 = ModifiedRate(rate,
                          new_reactants=["he4", "c12"],
                          new_products=["o16"], update_screening=False)

        assert (r1.eval(2e9, rho=1e8, comp=comp, screen_func=chugunov_2007) ==
                rate.eval(2e9, rho=1e8, comp=comp, screen_func=chugunov_2007))

        # Test eval against the original rate with updating ion_screen,
        # i.e. uses the new_reactants for screening
        r2 = ModifiedRate(rate,
                          new_reactants=["he4", "c12"],
                          new_products=["o16"], update_screening=True)

        # Sync ion_screen and screening_pair between the ModifiedRate
        # and the original rate then the results should match.
        rate.ion_screen = r2.ion_screen
        rate.screening_pairs = r2.screening_pairs

        assert (r2.eval(2e9, rho=1e8, comp=comp, screen_func=chugunov_2007) ==
                rate.eval(2e9, rho=1e8, comp=comp, screen_func=chugunov_2007))
