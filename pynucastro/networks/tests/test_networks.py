# unit tests for rates
import pynucastro.networks as networks
import pynucastro.rates as rates

from pytest import approx


class TestComposition(object):
    @classmethod
    def setup_class(cls):
        """ this is run once for each class before any tests """
        pass

    @classmethod
    def teardown_class(cls):
        """ this is run once for each class after all tests """
        pass

    def setup_method(self):
        """ this is run before each test """
        self.nuclei = [rates.Nucleus("h1"),
                       rates.Nucleus("he4"),
                       rates.Nucleus("c12"),
                       rates.Nucleus("o16"),
                       rates.Nucleus("n14"),
                       rates.Nucleus("ca40")]

        self.comp = networks.Composition(self.nuclei)

    def teardown_method(self):
        """ this is run after each test """
        self.tf = None

    def test_solar(self):
        self.comp.set_solar_like()

        sum = 0.0
        for k in self.comp.X:
            sum += self.comp.X[k]

        assert sum == approx(1.0)
        assert self.comp.X[rates.Nucleus("h1")] == approx(0.7)

    def test_set_all(self):
        val = 1.0/len(self.nuclei)
        self.comp.set_all(1.0/len(self.nuclei))
        for n in self.nuclei:
            assert self.comp.X[n] == val

    def test_set_nuc(self):
        n = self.nuclei[0]
        self.comp.set_nuc(n.raw, 0.55)
        assert self.comp.X[n] == 0.55

    def test_get_molar(self):
        self.comp.set_solar_like(Z=0.02)
        molar = self.comp.get_molar()
        assert molar[rates.Nucleus("he4")] == approx((0.3-0.02)/4.0)


class TestRateCollection(object):
    @classmethod
    def setup_class(cls):
        """ this is run once for each class before any tests """
        pass

    @classmethod
    def teardown_class(cls):
        """ this is run once for each class after all tests """
        pass

    def setup_method(self):
        """ this is run before each test """
        files = ["c12-pg-n13-ls09",
                 "c13-pg-n14-nacr",
                 "n13--c13-wc12",
                 "n13-pg-o14-lg06",
                 "n14-pg-o15-im05",
                 "n15-pa-c12-nacr",
                 "o14--n14-wc12",
                 "o15--n15-wc12"]
        self.rc = networks.RateCollection(files)

        self.p = rates.Nucleus("p")
        self.he4 = rates.Nucleus("he4")
        self.c12 = rates.Nucleus("c12")
        self.c13 = rates.Nucleus("c13")
        self.n13 = rates.Nucleus("n13")
        self.n14 = rates.Nucleus("n14")
        self.n15 = rates.Nucleus("n15")
        self.o14 = rates.Nucleus("o14")
        self.o15 = rates.Nucleus("o15")

    def teardown_method(self):
        """ this is run after each test """
        self.tf = None

    def test_nuclei(self):
        nuc = self.rc.get_nuclei()
        assert nuc == [self.p, self.he4, self.c12, self.c13,
                       self.n13, self.n14, self.n15, self.o14, self.o15]

    def test_eval(self):
        c = networks.Composition(self.rc.unique_nuclei)
        c.set_solar_like()

        rates = {"c12 + p --> n13": 4.3825344233265815e-05,
                 "c13 + p --> n14": 0.00012943869407433355,
                 "n13 --> c13": 2.5475016632596765e-07,
                 "n13 + p --> o14": 4.8517620910445875e-06,
                 "n14 + p --> o15": 9.8137074572314962e-07,
                 "n15 + p --> he4 + c12": 0.087518552257659241,
                 "o14 --> n14": 2.0036691481625654e-06,
                 "o15 --> n15": 1.0822012944765837e-06}

        rv = self.rc.evaluate_rates(1.e4, 1.e8, c)

        for r in rv:
            assert rv[r] == approx(rates[str(r)])

    def test_overview(self):

        ostr = """
  p
    consumed by:
       c12 + p --> n13
       c13 + p --> n14
       n13 + p --> o14
       n14 + p --> o15
       n15 + p --> he4 + c12
    produced by:

  he4
    consumed by:
    produced by:
       n15 + p --> he4 + c12

  c12
    consumed by:
       c12 + p --> n13
    produced by:
       n15 + p --> he4 + c12

  c13
    consumed by:
       c13 + p --> n14
    produced by:
       n13 --> c13

  n13
    consumed by:
       n13 --> c13
       n13 + p --> o14
    produced by:
       c12 + p --> n13

  n14
    consumed by:
       n14 + p --> o15
    produced by:
       c13 + p --> n14
       o14 --> n14

  n15
    consumed by:
       n15 + p --> he4 + c12
    produced by:
       o15 --> n15

  o14
    consumed by:
       o14 --> n14
    produced by:
       n13 + p --> o14

  o15
    consumed by:
       o15 --> n15
    produced by:
       n14 + p --> o15
"""
        assert self.rc.network_overview().replace(" ","").strip() == ostr.replace(" ","").strip()
