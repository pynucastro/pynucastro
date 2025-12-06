"""This is a simple test on the different IDs we have for rates.  It
serves as a comparison so we can remember what the different IDs mean.

"""

import pytest

import pynucastro as pyna
from pynucastro.rates.alternate_rates import DeBoerC12agO16


class TestRateIds:
    @pytest.fixture(scope="class")
    def approx_rate(self, reaclib_library):
        mylib = reaclib_library.linking_nuclei(["mg24", "al27", "si28",
                                                "p31", "s32", "he4", "p"])
        pynet = pyna.PythonNetwork(libraries=[mylib], verbose=True)
        pynet.make_ap_pg_approx()
        pynet.remove_nuclei(["al27", "p31"])
        ra = pynet.get_rates()[0]
        return ra

    def test_reaclib_rate(self, reaclib_library):
        r = reaclib_library.get_rate_by_name("c12(a,g)o16")

        assert r.rid == "C12 + He4 --> O16"
        assert r.id == "C12 + He4 --> O16 <nac2_reaclib__>"
        assert r.fname == "He4_C12__O16"
        assert r.cname() == "He4_C12_to_O16"

    def test_alternate_rate(self):
        r = DeBoerC12agO16()

        assert r.rid == "C12 + He4 --> O16"
        assert r.id == "C12 + He4 --> O16 <debo_reaclib__>"
        assert r.fname == "C12_He4__O16"
        assert r.cname() == "C12_He4_to_O16"

    def test_tabular_rate(self, langanke_library, ffn_library):

        # get the same rare from 2 different libraries.  They
        # should really be distinguishable.

        r1 = langanke_library.get_rate_by_name("ni56(,)co56")

        assert r1.rid == "Ni56 --> Co56"
        assert r1.id == "Ni56 --> Co56 <tabular_tabular>"
        assert r1.fname == "Ni56__Co56"
        assert r1.cname() == "Ni56_to_Co56"

        r2 = ffn_library.get_rate_by_name("ni56(,)co56")

        assert r2.rid == "Ni56 --> Co56"
        assert r2.id == "Ni56 --> Co56 <tabular_tabular>"
        assert r2.fname == "Ni56__Co56"
        assert r2.cname() == "Ni56_to_Co56"

        assert r1 == r2

    def test_derived_rate(self, reaclib_library):
        r = reaclib_library.get_rate_by_name("c12(a,g)o16")

        rr = pyna.DerivedRate(r, compute_Q=True, use_pf=True)

        assert rr.rid == "O16 --> He4 + C12"
        assert rr.id == "O16 --> He4 + C12 <derived_reaclib__derived_from_inverse>"
        assert rr.fname == "O16__He4_C12__derived"
        assert rr.cname() == "O16_to_He4_C12_derived"

    def test_approximate_rate(self, approx_rate):

        ra = approx_rate

        assert ra.rid == "Mg24 + He4 --> Si28"
        assert ra.id == "Mg24 + He4 --> Si28 <approx>"
        assert ra.fname == "Mg24_He4__Si28__approx"
        assert ra.cname() == "Mg24_He4_to_Si28_approx"
