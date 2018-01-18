# unit tests for rates
import math

import pyreaclib.networks as networks
import pyreaclib.rates as rates

from pytest import approx

class TestFortranNetwork(object):
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
                 "c12-c12p-na23-cf88",
                 "c13-pg-n14-nacr",
                 "n13--c13-wc12",
                 "n13-pg-o14-lg06",
                 "n14-pg-o15-im05",
                 "n15-pa-c12-nacr",
                 "o14--n14-wc12",
                 "o15--n15-wc12"]
        self.fn = networks.StarKillerNetwork(files)

        self.p = rates.Nucleus("p")
        self.he4 = rates.Nucleus("he4")
        self.c12 = rates.Nucleus("c12")
        self.c13 = rates.Nucleus("c13")
        self.n13 = rates.Nucleus("n13")
        self.n14 = rates.Nucleus("n14")
        self.n15 = rates.Nucleus("n15")
        self.o14 = rates.Nucleus("o14")
        self.o15 = rates.Nucleus("o15")
        self.na23 = rates.Nucleus("na23")

    def teardown_method(self):
        """ this is run after each test """
        self.tf = None


    def test_ydot_string(self):
        """ test the ydot string for this network """

        # Test a reaction of the type A + B -> C where A != B
        answer = "dens * Y(jp) * Y(jc12) * screened_rates(i_scor, k_c12_pg_n13) * screened_rates(i_rate, k_c12_pg_n13)"

        ydot = self.fn.ydot_string(self.fn.rates[0])
        assert ydot == answer

        # Test a reaction of the type A + B -> C + D where A == B
        answer = "5.00000000000000d-01 * dens * Y(jc12)**2 * screened_rates(i_scor, k_c12_c12p_na23) * screened_rates(i_rate, k_c12_c12p_na23)"

        ydot = self.fn.ydot_string(self.fn.rates[1])
        assert ydot == answer

class TestReaclibChapterNetwork(object):
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

        files = ["b17-nnn-c14-wc12",
                 "he3-he3pp-he4-nacr",
                 "he4-aag-c12-fy05",
                 "he4-npahe3-li7-mafo",
                 "he4-pphe3-he3-nacr",
                 "he6-gnn-he4-cf88",
                 "li7-tnna-he4-mafo",
                 "n--p-wc12",
                 "p-ng-d-an06",
                 "t-gn-d-nk06",
                 "t-pn-he3-de04"]

        self.fn = networks.StarKillerNetwork(files)

        self.n = rates.Nucleus("n")
        self.p = rates.Nucleus("p")
        self.d = rates.Nucleus("d")
        self.t = rates.Nucleus("t")
        self.he3 = rates.Nucleus("he3")
        self.he4 = rates.Nucleus("he4")
        self.he6 = rates.Nucleus("he6")
        self.li7 = rates.Nucleus("li7")
        self.b17 = rates.Nucleus("b17")
        self.c12 = rates.Nucleus("c12")
        self.c14 = rates.Nucleus("c14")

    def teardown_method(self):
        """ this is run after each test """
        self.tf = None

    def test_nuclei(self):
        """ test the nuclei are correctly identified for this network """
        assert self.n.N == 1
        assert self.n.Z == 0
        assert self.n.A == 1

        assert self.p.N == 0
        assert self.p.Z == 1
        assert self.p.A == 1

        assert self.d.N == 1
        assert self.d.Z == 1
        assert self.d.A == 2

        assert self.t.N == 2
        assert self.t.Z == 1
        assert self.t.A == 3

        assert self.he3.N == 1
        assert self.he3.Z == 2
        assert self.he3.A == 3

        assert self.he4.N == 2
        assert self.he4.Z == 2
        assert self.he4.A == 4

        assert self.he6.N == 4
        assert self.he6.Z == 2
        assert self.he6.A == 6

        assert self.li7.N == 4
        assert self.li7.Z == 3
        assert self.li7.A == 7

        assert self.b17.N == 12
        assert self.b17.Z == 5
        assert self.b17.A == 17

        assert self.c12.N == 6
        assert self.c12.Z == 6
        assert self.c12.A == 12

        assert self.c14.N == 8
        assert self.c14.Z == 6
        assert self.c14.A == 14
