# unit tests for screening with approximate rates
# we'll test both symmetric and "normal" screening

import pynucastro as pyna


class TestApproxScreening:
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

        rl = pyna.ReacLibLibrary()
        mynet = rl.linking_nuclei(["p", "he4", "mg24",
                                   "al27", "si28", "p31", "s32"])

        pynet_symmetric = pyna.PythonNetwork(libraries=[mynet], symmetric_screening=True)
        pynet_symmetric.make_ap_pg_approx()
        pynet_symmetric.remove_nuclei(["al27", "p31"])
        self.pynet_symmetric = pynet_symmetric

        pynet = pyna.PythonNetwork(libraries=[mynet])
        pynet.make_ap_pg_approx()
        pynet.remove_nuclei(["al27", "p31"])
        self.pynet = pynet

    def teardown_method(self):
        """ this is run after each test """
        self.pynet = None
        self.pynet_symmetric = None

    def test_symmetric_screening(self):
        screening_map = self.pynet_symmetric.get_screening_map()

        # all of the reaclib rates that are 2 body (and not "n") should have
        # be in the screening map

        for r in self.pynet_symmetric.reaclib_rates:
            if r.reverse:
                nucs = [q for q in r.products if q.Z != 0]
                if len(nucs) == 1:
                    continue
                assert nucs[0] in r.symmetric_screen and nucs[1] in r.symmetric_screen
                r_scn = [q for q in screening_map if q.n1 in r.symmetric_screen and q.n2 in r.symmetric_screen]
                assert len(r_scn) == 1
            else:
                nucs = [q for q in r.reactants if q.Z != 0]
                if len(nucs) == 1:
                    continue
                assert nucs[0] in r.symmetric_screen and nucs[1] in r.symmetric_screen
                r_scn = [q for q in screening_map if q.n1 in r.symmetric_screen and q.n2 in r.symmetric_screen]
                assert len(r_scn) == 1

    def test_screening(self):
        screening_map = self.pynet.get_screening_map()

        # all of the reaclib rates that are 2 body (and not "n") should have
        # be in the screening map

        for r in self.pynet.reaclib_rates:
            nucs = [q for q in r.reactants if q.Z != 0]
            if len(nucs) == 1:
                continue
            assert nucs[0] in r.ion_screen and nucs[1] in r.ion_screen
            r_scn = [q for q in screening_map if q.n1 in r.ion_screen and q.n2 in r.ion_screen]
            assert len(r_scn) == 1

