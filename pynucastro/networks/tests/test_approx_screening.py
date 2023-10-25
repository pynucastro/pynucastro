# unit tests for screening with approximate rates
# we'll test both symmetric and "normal" screening

import pytest

import pynucastro as pyna
from pynucastro.screening import get_screening_map


class TestApproxScreening:
    @pytest.fixture(scope="class")
    def mynet(self, reaclib_library):
        return reaclib_library.linking_nuclei(["p", "he4", "mg24",
                                               "al27", "si28", "p31", "s32"])

    @pytest.fixture(scope="class")
    def pynet_symmetric(self, mynet):
        pynet_symmetric = pyna.PythonNetwork(libraries=[mynet], symmetric_screening=True)
        pynet_symmetric.make_ap_pg_approx()
        pynet_symmetric.remove_nuclei(["al27", "p31"])
        return pynet_symmetric

    @pytest.fixture(scope="class")
    def pynet(self, mynet):
        pynet = pyna.PythonNetwork(libraries=[mynet])
        pynet.make_ap_pg_approx()
        pynet.remove_nuclei(["al27", "p31"])
        return pynet

    def test_symmetric_screening(self, pynet_symmetric):
        screening_map = get_screening_map(pynet_symmetric.get_rates(),
                                          symmetric_screening=pynet_symmetric.symmetric_screening)

        # all of the reaclib rates that are 2 body (and not "n") should have
        # be in the screening map

        for r in pynet_symmetric.reaclib_rates:
            if r.Q < 0:
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

    def test_screening(self, pynet):
        screening_map = get_screening_map(pynet.get_rates(),
                                          symmetric_screening=pynet.symmetric_screening)

        # all of the reaclib rates that are 2 body (and not "n") should have
        # be in the screening map

        for r in pynet.reaclib_rates:
            nucs = [q for q in r.reactants if q.Z != 0]
            if len(nucs) == 1:
                continue
            assert nucs[0] in r.ion_screen and nucs[1] in r.ion_screen
            r_scn = [q for q in screening_map if q.n1 in r.ion_screen and q.n2 in r.ion_screen]
            assert len(r_scn) == 1
