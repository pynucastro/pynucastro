"""This is a simple test on the different IDs we have for rates.  It
serves as a comparison so we can remember what the different IDs mean.

"""

import pytest

import pynucastro as pyna
from pynucastro.rates.alternate_rates import DeBoerC12agO16, IliadisO16pgF17
from pynucastro.rates.known_duplicates import ALLOWED_DUPLICATES


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
        r1 = reaclib_library.get_rate_by_name("c12(a,g)o16")

        assert r1.rid == "C12 + He4 --> O16"
        assert r1.id == "C12 + He4 --> O16 <reaclib_nac2>"
        assert r1.fname == "He4_C12_to_O16_reaclib"

    def test_alternate_rate(self):
        r1 = DeBoerC12agO16()

        assert r1.rid == "C12 + He4 --> O16"
        assert r1.id == "C12 + He4 --> O16 <deboer_debo>"
        assert r1.fname == "C12_He4_to_O16_deboer"

        r2 = IliadisO16pgF17()

        assert r2.rid == "O16 + p --> F17"
        assert r2.id == "O16 + p --> F17 <iliadis_iliadis2022>"
        assert r2.fname == "O16_p_to_F17_iliadis"

    def test_tabular_rate(self, langanke_library, ffn_library):

        # get the same rare from 2 different libraries.  They
        # should really be distinguishable.

        r1 = langanke_library.get_rate_by_name("ni56(,)co56")

        assert r1.rid == "Ni56 --> Co56"
        assert r1.id == "Ni56 --> Co56 <weaktab_langanke>"
        assert r1.fname == "Ni56_to_Co56_weaktab"

        r2 = ffn_library.get_rate_by_name("ni56(,)co56")

        assert r2.rid == "Ni56 --> Co56"
        assert r2.id == "Ni56 --> Co56 <weaktab_ffn>"
        assert r2.fname == "Ni56_to_Co56_weaktab"

        assert r1 == r2

    def test_derived_rate(self, reaclib_library):
        r = reaclib_library.get_rate_by_name("c12(a,g)o16")

        rr = pyna.DerivedRate(r, compute_Q=True, use_pf=True)

        assert rr.rid == "O16 --> He4 + C12"
        assert rr.id == "O16 --> He4 + C12 <derived_nac2>"
        assert rr.fname == "O16_to_He4_C12_derived"

    def test_approximate_rate(self, approx_rate):

        ra = approx_rate

        assert ra.rid == "Mg24 + He4 --> Si28"
        assert ra.id == "Mg24 + He4 --> Si28 <approx>"
        assert ra.fname == "Mg24_He4_to_Si28_approx"

        cr = ra.get_child_rates()[1]

        assert cr.rid == "Mg24 + He4 --> p + Al27"
        assert cr.id == "Mg24 + He4 --> p + Al27 <removed_il10>"
        assert cr.fname == "He4_Mg24_to_p_Al27_removed"

    def test_duplicate_rates(self, full_library):
        # Make sure duplicate sets in ALLOWED_DUPLICATES have different fname and id

        for dupe_set in ALLOWED_DUPLICATES:
            dupe_rates = []
            for dupe_entry in dupe_set:
                rate_id = dupe_entry.split(":", 1)[1].strip()
                dupe_rates.append(full_library.get_rate(rate_id))

            fnames = [r.fname for r in dupe_rates]
            rate_ids = [r.id for r in dupe_rates]

            # Make sure fnames and rate ids are unique.
            assert len(fnames) == len(set(fnames))
            assert len(rate_ids) == len(set(rate_ids))
