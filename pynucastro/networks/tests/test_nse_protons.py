import pynucastro as pyna


class TestNseProtons:

    def test_nse_protons(self, capsys):

        nuc = ["p", "he4", "ti44", "v47", "cr48", "mn51",
               "fe52", "co55", "ni56"]

        # note: we can't use the reaclib_library fixture we've
        # setup in pytest because we will be modifying rates,
        # which will be reflected in the fixture and affect
        # later tests
        rl = pyna.ReacLibLibrary()

        lib = rl.linking_nuclei(nuc)

        rc = pyna.RateCollection(libraries=lib)

        rc.make_nse_protons(51)

        captured = capsys.readouterr()

        output = """modifying p_Mn51__Fe52 to use NSE protons
modifying Fe52__p_Mn51 to use NSE protons
modifying p_Co55__Ni56 to use NSE protons
modifying Ni56__p_Co55 to use NSE protons
modifying He4_Fe52__p_Co55 to use NSE protons
modifying p_Co55__He4_Fe52 to use NSE protons
"""

        assert captured.out == output
