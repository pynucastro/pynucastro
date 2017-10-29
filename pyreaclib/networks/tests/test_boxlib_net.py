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
                 "c13-pg-n14-nacr",
                 "n13--c13-wc12",
                 "n13-pg-o14-lg06",
                 "n14-pg-o15-im05",
                 "n15-pa-c12-nacr",
                 "o14--n14-wc12",
                 "o15--n15-wc12"]
        self.fn = networks.BoxLibNetwork(files)

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


    def test_ydot_string(self):
        """ test the ydot string for this network """

        answer_a = "dens * Y(jc12) * Y(jp) * screened_rates(i_scor, k_c12_pg_n13) * screened_rates(i_rate, k_c12_pg_n13)"

        answer_b = "dens * Y(jp) * Y(jc12) * screened_rates(i_scor, k_c12_pg_n13) * screened_rates(i_rate, k_c12_pg_n13)"

        ydot = self.fn.ydot_string(self.fn.rates[0])
        assert ydot == answer_a or ydot == answer_b


