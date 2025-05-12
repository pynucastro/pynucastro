# unit tests for rates
import math

import pytest
from pytest import approx

from pynucastro import Composition, Rate, rates
from pynucastro.nucdata import Nucleus
from pynucastro.rates import BaryonConservationError


class TestTfactors:
    @pytest.fixture(scope="class")
    def tf(self):
        return rates.Tfactors(2.e9)

    def test_tfactors(self, tf):
        assert tf.T9 == approx(2.0)
        assert tf.T9i == approx(0.5)
        assert tf.T913i == approx(0.5**(1./3.))
        assert tf.T913 == approx(2.0**(1./3.))
        assert tf.T953 == approx(2.0**(5./3.))
        assert tf.lnT9 == approx(math.log(2.0))


class TestRate:

    @classmethod
    def setup_class(cls):
        """ this is run once for each class before any tests """

    @classmethod
    def teardown_class(cls):
        """ this is run once for each class after all tests """

    def setup_method(self):
        """ this is run before each test """

        # chapter-1
        self.rate1 = rates.load_rate("o15--n15-wc12")

        # chapter-2
        self.rate2 = rates.load_rate("t-gn-d-nk06")

        # chapter-3
        self.rate3 = rates.load_rate("he6-gnn-he4-cf88")

        # chapter-4
        self.rate4 = rates.load_rate("c12-ag-o16-nac2")

        # chapter-5
        self.rate5 = rates.load_rate("n15-pa-c12-nacr")

        # chapter-6
        self.rate6 = rates.load_rate("he3-he3pp-he4-nacr")

        # chapter-7
        self.rate7 = rates.load_rate("li7-tnna-he4-mafo")

        # chapter-8
        self.rate8 = rates.load_rate("he4-aag-c12-fy05")

        # chapter-8, historical format (same rate as chapter-9)
        self.rate8_hist = rates.load_rate("he4-pphe3-he3-nacr-historical")

        # chapter-9
        self.rate9 = rates.load_rate("he4-pphe3-he3-nacr")

        # chapter-10
        self.rate10 = rates.load_rate("he4-npahe3-li7-mafo")

        # chapter-11
        self.rate11 = rates.load_rate("b17-nnn-c14-wc12")

        self.n = Nucleus("n")

        self.p = Nucleus("p")
        self.h1 = Nucleus("H1")
        self.d = Nucleus("d")
        self.h3 = Nucleus("H3")

        self.he3 = Nucleus("He3")
        self.he4 = Nucleus("He4")
        self.he6 = Nucleus("He6")

        self.li7 = Nucleus("Li7")

        self.b17 = Nucleus("B17")

        self.c12 = Nucleus("C12")
        self.c14 = Nucleus("C14")

        self.n15 = Nucleus("N15")

        self.o15 = Nucleus("O15")
        self.o16 = Nucleus("O16")

        self.ni56 = Nucleus("Ni56")
        self.u238 = Nucleus("U238")
        self.he4_also = Nucleus("he4")

    def teardown_method(self):
        """ this is run after each test """

    def test_source(self):
        assert self.rate1.source["Year"] == "2012"

    def test_reactants(self):

        # o15--n15-wc12
        assert self.rate1.reactants[0] == self.o15
        assert len(self.rate1.reactants) == 1

        # t-gn-d-nk06
        assert self.rate2.reactants[0] == self.h3
        assert len(self.rate2.reactants) == 1

        # he6-gnn-he4-cf88
        assert self.rate3.reactants[0] == self.he6
        assert len(self.rate3.reactants) == 1

        # c12-ag-o16-nac2
        assert self.rate4.reactants[0] == self.he4
        assert self.rate4.reactants[1] == self.c12
        assert len(self.rate4.reactants) == 2

        # n15-pa-c12-nacr
        assert self.rate5.reactants[0] == self.h1
        assert self.rate5.reactants[1] == self.n15
        assert len(self.rate5.reactants) == 2

        # he3-he3pp-he4-nacr
        assert self.rate6.reactants[0] == self.he3
        assert self.rate6.reactants[1] == self.he3
        assert len(self.rate6.reactants) == 2

        # li7-tnna-he4-mafo
        assert self.rate7.reactants[0] == self.h3
        assert self.rate7.reactants[1] == self.li7
        assert len(self.rate7.reactants) == 2

        # he4-aag-c12-fy05
        assert self.rate8.reactants[0] == self.he4
        assert self.rate8.reactants[1] == self.he4
        assert self.rate8.reactants[2] == self.he4
        assert len(self.rate8.reactants) == 3

        # he4-pphe3-he3-nacr-historical
        assert self.rate8_hist.reactants[0] == self.p
        assert self.rate8_hist.reactants[1] == self.h1
        assert self.rate8_hist.reactants[2] == self.he4
        assert len(self.rate8_hist.reactants) == 3

        # he4-pphe3-he3-nacr
        assert self.rate9.reactants[0] == self.p
        assert self.rate9.reactants[1] == self.h1
        assert self.rate9.reactants[2] == self.he4
        assert len(self.rate9.reactants) == 3

        # he4-npahe3-li7-mafo
        assert self.rate10.reactants[0] == self.n
        assert self.rate10.reactants[1] == self.h1
        assert self.rate10.reactants[2] == self.he4
        assert self.rate10.reactants[3] == self.he4
        assert len(self.rate10.reactants) == 4

        # b17-nnn-c14-wc12
        assert self.rate11.reactants[0] == self.b17
        assert len(self.rate11.reactants) == 1

    def test_products(self):
        assert self.rate4.products[0] == self.o16
        assert self.rate8.products[0] == self.c12
        assert len(self.rate8.products) == 1

        # he4-pphe3-he3-nacr-historical
        assert self.rate8_hist.products[0] == self.he3
        assert self.rate8_hist.products[1] == self.he3
        assert len(self.rate8_hist.products) == 2

        # he4-pphe3-he3-nacr
        assert self.rate9.products[0] == self.he3
        assert self.rate9.products[1] == self.he3
        assert len(self.rate9.products) == 2

    def test_count(self):
        assert self.rate3.reactant_count(Nucleus("he6")) == 1
        assert self.rate3.product_count(Nucleus("he6")) == 0

        assert self.rate5.reactant_count(Nucleus("n15")) == 1
        assert self.rate5.reactant_count(Nucleus("p")) == 1
        assert self.rate5.reactant_count(Nucleus("a")) == 0
        assert self.rate5.reactant_count(Nucleus("c12")) == 0

        assert self.rate5.product_count(Nucleus("n15")) == 0
        assert self.rate5.product_count(Nucleus("p")) == 0
        assert self.rate5.product_count(Nucleus("a")) == 1
        assert self.rate5.product_count(Nucleus("c12")) == 1

        assert self.rate8.reactant_count(Nucleus("a")) == 3

    def test_prefactor(self):
        assert self.rate4.prefactor == 1.0
        assert self.rate8.prefactor == approx(0.16666666)

    def test_rate_exponent(self):
        assert self.rate8.get_rate_exponent(1.e8) == approx(40.9106396)

    def test_eval(self):
        assert self.rate8.eval(1.e8) == approx(2.0403192412842946e-24, rel=1.e-6, abs=1.e-40)

    def test_eval_deriv(self):
        T0 = 1.e8
        eps = 1.e-8

        # compare finite diff to analytic diff

        # rate4
        diff = (self.rate4.eval(T0*(1.0+eps)) - self.rate4.eval(T0)) / (T0 * eps)
        err = abs(diff - self.rate4.eval_deriv(T0)) / diff

        assert err < 1.e-6

        # rate5
        diff = (self.rate5.eval(T0*(1.0+eps)) - self.rate5.eval(T0)) / (T0 * eps)
        err = abs(diff - self.rate5.eval_deriv(T0)) / diff

        assert err < 1.e-6

        # rate6
        diff = (self.rate6.eval(T0*(1.0+eps)) - self.rate6.eval(T0)) / (T0 * eps)
        err = abs(diff - self.rate6.eval_deriv(T0)) / diff

        assert err < 1.e-6

    def test_comparison(self):
        assert self.rate1 > self.rate2
        assert self.rate1 > self.rate4
        assert self.rate8 > self.rate9

    def test_weak(self):
        assert self.rate1.weak
        assert not self.rate2.weak

    def test_screen(self):
        assert not self.rate1.ion_screen
        assert self.rate4.ion_screen == [Nucleus("he4"), Nucleus("c12")]
        assert self.rate8.ion_screen == 3*[Nucleus("he4")]

    def test_heaviest_lightest(self):
        assert self.rate4.heaviest() == Nucleus("o16")
        assert self.rate4.lightest() == Nucleus("he4")
        assert self.rate2.lightest() == Nucleus("n")
        assert self.rate2.heaviest() == Nucleus("t")

    def test_identical_particle_factor(self):
        assert self.rate8.prefactor == approx(0.16666667)

        self.rate8.use_identical_particle_factor = False
        self.rate8._set_rhs_properties()  # pylint: disable=protected-access

        assert self.rate8.prefactor == 1.0

    def test_stoichiometry(self, reaclib_library):
        assert repr(self.rate4) == "C12 + He4 ⟶ O16 + 𝛾"

        # create a separate version since rates are mutable
        c12ag = reaclib_library.get_rate_by_name("c12(a,g)o16")
        c12ag.stoichiometry = {Nucleus("he4"): 1.5,
                               Nucleus("c12"): 1,
                               Nucleus("o16"): 1}
        c12ag._set_print_representation()  # pylint: disable=protected-access

        assert repr(c12ag) == "C12 + 1.5 He4 ⟶ O16"
        assert c12ag.rid == "C12 + 1.5 He4 --> O16"

        assert c12ag.reactant_count(Nucleus("he4")) == 1.5

        # restore it so the library is unchanged
        c12ag.stoichiometry = None

    def test_stoichiometry_3alpha(self, reaclib_library):

        three_alpha = reaclib_library.get_rate_by_name("he4(aa,g)c12")
        assert repr(three_alpha) == "3 He4 ⟶ C12 + 𝛾"
        assert three_alpha.rid == "3 He4 --> C12"

        three_alpha.stoichiometry = {Nucleus("he4"): 4,
                                     Nucleus("c12"): 1}
        three_alpha._set_print_representation()  # pylint: disable=protected-access

        assert repr(three_alpha) == "4 He4 ⟶ C12"
        assert three_alpha.rid == "4 He4 --> C12"

        assert three_alpha.reactant_count(Nucleus("he4")) == 4

        three_alpha.stoichiometry = None

    def test_stoichiometry_dict(self, reaclib_library):

        c12c12 = reaclib_library.get_rate_by_name("c12(c12,a)ne20")

        c12c12.stoichiometry = {Nucleus("he4"): 4}
        c12c12._set_print_representation()  # pylint: disable=protected-access

        assert repr(c12c12) == "C12 + C12 ⟶ 4 He4 + Ne20"
        assert c12c12.rid == "C12 + C12 --> 4 He4 + Ne20"

        assert c12c12.reactant_count(Nucleus("he4")) == 0
        assert c12c12.product_count(Nucleus("he4")) == 4

        assert c12c12.reactant_count(Nucleus("c12")) == 2

        assert c12c12.product_count(Nucleus("ne20")) == 1

        c12c12.stoichiometry = None

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

    def a_a_ag_c12(self, reaclib_library):
        """
        Here we test the inverse rate, computed by the use of detailed balance
        of a:

        A + B -> C

        reaction type.
        """

        a_a_ag_c12 = reaclib_library.get_rate_by_name("he4(aa,g)c12")
        c12_ga_a_a_reaclib = reaclib_library.get_rate_by_name("c12(g,aa)he4")
        c12_ga_a_a_derived = rates.DerivedRate(rate=a_a_ag_c12, compute_Q=False, use_pf=False)

        assert c12_ga_a_a_reaclib.eval(T=2.0e9) == approx(c12_ga_a_a_derived.eval(T=2.0e9), rel=2e-4)

    def test_a_a_ag_c12_with_pf(self, reaclib_library):
        """
        This function test the correct rate value if we take in consideration the partition
        functions on the range 1.0e9 to 100.0e9
        """

        a_a_ag_c12 = reaclib_library.get_rate_by_name("he4(aa,g)c12")
        c12_ga_a_a_derived = rates.DerivedRate(rate=a_a_ag_c12, compute_Q=False, use_pf=True)

        with pytest.warns(UserWarning, match="C12 partition function is not supported by tables"):
            rval = c12_ga_a_a_derived.eval(T=2.0e9)
        assert rval == approx(2.8953989705969484e-07)

    def test_a_a_ag_c12_with_Q(self, reaclib_library):
        """
        This function test the correct rate value if we take in consideration the
        exact values of atomic nuclear weight in order to compute the Q capture value
        of the reaction rate.
        """

        a_a_ag_c12 = reaclib_library.get_rate_by_name("he4(aa,g)c12")
        c12_ga_a_a_derived = rates.DerivedRate(rate=a_a_ag_c12, compute_Q=True, use_pf=False)

        assert c12_ga_a_a_derived.eval(T=2.0e9) == approx(2.899642192191721e-07)


class TestWeakRates:
    @pytest.fixture(scope="class")
    def rate1(self):
        return rates.TabularRate("suzuki-o18--f18-toki")

    @pytest.fixture(scope="class")
    def rate2(self):
        return rates.TabularRate("suzuki-na22--ne22-toki")

    @pytest.fixture(scope="class")
    def rate3(self):
        return rates.TabularRate("langanke-sc45--ca45-toki")

    @pytest.fixture(scope="class")
    def rate4(self):
        return rates.TabularRate("langanke-ti45--sc45-toki")

    @pytest.fixture(scope="class")
    def rate5(self):
        return rates.TabularRate("langanke-v45--ti45-toki")

    @pytest.fixture(scope="class")
    def rate6(self):
        return rates.TabularRate("langanke-ca45--sc45-toki")

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
    def rate(self):
        return rates.load_rate("c12-c12n-mg23-cf88")

    def test_modify(self, rate):

        rate.modify_products("mg24")

        assert rate.Q == approx(13.933578000000125)
        assert rate.products == [Nucleus("mg24")]
        assert rate.modified
