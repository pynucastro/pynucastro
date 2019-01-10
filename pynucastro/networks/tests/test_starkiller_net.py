# unit tests for rates
import pynucastro.networks as networks
import pynucastro.rates as rates

import io


class TestStarKillerNetwork(object):
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
        files = ["c12-c12a-ne20-cf88",
                 "c12-c12n-mg23-cf88",
                 "c12-c12p-na23-cf88",
                 "c12-ag-o16-nac2",
                 "na23--ne23-toki",
                 "ne23--na23-toki",
                 "n--p-wc12"]

        self.fn = networks.StarKillerNetwork(files)

    def teardown_method(self):
        """ this is run after each test """
        self.tf = None

    def cromulent_ftag(self, ftag, answer, n_indent=1):
        """ check to see if function ftag returns answer """

        output = io.StringIO()
        ftag(n_indent, output)
        result = output.getvalue() == answer
        output.close()
        return result

    def test_nrates(self):
        """ test the _nrates function """

        answer = '  integer, parameter :: nrates = 7\n'
        assert self.cromulent_ftag(self.fn._nrates, answer, n_indent=1)

    def test_nrat_reaclib(self):
        """ test the _nrat_reaclib function """

        answer = ('  integer, parameter :: nrat_reaclib = 5\n' +
                  '  integer, parameter :: number_reaclib_sets = 6\n')
        assert self.cromulent_ftag(self.fn._nrat_reaclib, answer, n_indent=1)

    def test_nrat_tabular(self):
        """ test the _nrat_tabular function """

        answer = '  integer, parameter :: nrat_tabular = 2\n'
        assert self.cromulent_ftag(self.fn._nrat_tabular, answer, n_indent=1)

    def test_nspec(self):
        """ test the _nspec function """

        answer = '  integer, parameter :: nspec = 9\n'
        assert self.cromulent_ftag(self.fn._nspec, answer, n_indent=1)

    def test_nspec_evolve(self):
        """ test the _nspec_evolve function """

        answer = '  integer, parameter :: nspec_evolve = 9\n'
        assert self.cromulent_ftag(self.fn._nspec_evolve, answer, n_indent=1)

    def test_nrxn(self):
        """ test the _nrxn function """

        answer = ('  integer, parameter :: k_c12_c12__he4_ne20   = 1\n' +
                  '  integer, parameter :: k_c12_c12__n_mg23   = 2\n' +
                  '  integer, parameter :: k_c12_c12__p_na23   = 3\n' +
                  '  integer, parameter :: k_he4_c12__o16   = 4\n' +
                  '  integer, parameter :: k_n__p__weak__wc12   = 5\n' +
                  '  integer, parameter :: k_na23__ne23   = 6\n' +
                  '  integer, parameter :: k_ne23__na23   = 7\n')
        assert self.cromulent_ftag(self.fn._nrxn, answer, n_indent=1)

    def test_jion(self):
        """ test the _jion function """

        answer = ('  integer, parameter :: jn   = 1\n' +
                  '  integer, parameter :: jp   = 2\n' +
                  '  integer, parameter :: jhe4   = 3\n' +
                  '  integer, parameter :: jc12   = 4\n' +
                  '  integer, parameter :: jo16   = 5\n' +
                  '  integer, parameter :: jne20   = 6\n' +
                  '  integer, parameter :: jne23   = 7\n' +
                  '  integer, parameter :: jna23   = 8\n' +
                  '  integer, parameter :: jmg23   = 9\n')
        assert self.cromulent_ftag(self.fn._jion, answer, n_indent=1)

    def test_spec_names(self):
        """ test the _spec_names function """

        answer = ('    spec_names(jn)   = "neutron"\n' +
                  '    spec_names(jp)   = "hydrogen-1"\n' +
                  '    spec_names(jhe4)   = "helium-4"\n' +
                  '    spec_names(jc12)   = "carbon-12"\n' +
                  '    spec_names(jo16)   = "oxygen-16"\n' +
                  '    spec_names(jne20)   = "neon-20"\n' +
                  '    spec_names(jne23)   = "neon-23"\n' +
                  '    spec_names(jna23)   = "sodium-23"\n' +
                  '    spec_names(jmg23)   = "magnesium-23"\n')
        assert self.cromulent_ftag(self.fn._spec_names, answer, n_indent=2)

    def test_short_spec_names(self):
        """ test the _short_spec_names function """

        answer = ('    short_spec_names(jn)   = "n"\n' +
                  '    short_spec_names(jp)   = "h1"\n' +
                  '    short_spec_names(jhe4)   = "he4"\n' +
                  '    short_spec_names(jc12)   = "c12"\n' +
                  '    short_spec_names(jo16)   = "o16"\n' +
                  '    short_spec_names(jne20)   = "ne20"\n' +
                  '    short_spec_names(jne23)   = "ne23"\n' +
                  '    short_spec_names(jna23)   = "na23"\n' +
                  '    short_spec_names(jmg23)   = "mg23"\n')
        assert self.cromulent_ftag(self.fn._short_spec_names, answer, n_indent=2)

    def test_ebind(self):
        """ test the _ebind function """

        answer = ('    ebind_per_nucleon(jn)   = 0.00000000000000d+00\n' +
                  '    ebind_per_nucleon(jp)   = 0.00000000000000d+00\n' +
                  '    ebind_per_nucleon(jhe4)   = 7.07391500000000d+00\n' +
                  '    ebind_per_nucleon(jc12)   = 7.68014400000000d+00\n' +
                  '    ebind_per_nucleon(jo16)   = 7.97620600000000d+00\n' +
                  '    ebind_per_nucleon(jne20)   = 8.03224000000000d+00\n' +
                  '    ebind_per_nucleon(jne23)   = 7.95525600000000d+00\n' +
                  '    ebind_per_nucleon(jna23)   = 8.11149300000000d+00\n' +
                  '    ebind_per_nucleon(jmg23)   = 7.90111500000000d+00\n')
        assert self.cromulent_ftag(self.fn._ebind, answer, n_indent=2)

    def test_aion(self):
        """ test the _aion function """

        answer = ('    aion(jn)   = 1.00000000000000d+00\n' +
                  '    aion(jp)   = 1.00000000000000d+00\n' +
                  '    aion(jhe4)   = 4.00000000000000d+00\n' +
                  '    aion(jc12)   = 1.20000000000000d+01\n' +
                  '    aion(jo16)   = 1.60000000000000d+01\n' +
                  '    aion(jne20)   = 2.00000000000000d+01\n' +
                  '    aion(jne23)   = 2.30000000000000d+01\n' +
                  '    aion(jna23)   = 2.30000000000000d+01\n' +
                  '    aion(jmg23)   = 2.30000000000000d+01\n')
        assert self.cromulent_ftag(self.fn._aion, answer, n_indent=2)

    def test_zion(self):
        """ test the zion function """

        answer = ('    zion(jn)   = 0.00000000000000d+00\n' +
                  '    zion(jp)   = 1.00000000000000d+00\n' +
                  '    zion(jhe4)   = 2.00000000000000d+00\n' +
                  '    zion(jc12)   = 6.00000000000000d+00\n' +
                  '    zion(jo16)   = 8.00000000000000d+00\n' +
                  '    zion(jne20)   = 1.00000000000000d+01\n' +
                  '    zion(jne23)   = 1.00000000000000d+01\n' +
                  '    zion(jna23)   = 1.10000000000000d+01\n' +
                  '    zion(jmg23)   = 1.20000000000000d+01\n')
        assert self.cromulent_ftag(self.fn._zion, answer, n_indent=2)

    def test_nion(self):
        """ test the _nion function """

        answer = ('    nion(jn)   = 1.00000000000000d+00\n' +
                  '    nion(jp)   = 0.00000000000000d+00\n' +
                  '    nion(jhe4)   = 2.00000000000000d+00\n' +
                  '    nion(jc12)   = 6.00000000000000d+00\n' +
                  '    nion(jo16)   = 8.00000000000000d+00\n' +
                  '    nion(jne20)   = 1.00000000000000d+01\n' +
                  '    nion(jne23)   = 1.30000000000000d+01\n' +
                  '    nion(jna23)   = 1.20000000000000d+01\n' +
                  '    nion(jmg23)   = 1.10000000000000d+01\n')
        assert self.cromulent_ftag(self.fn._nion, answer, n_indent=2)


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
