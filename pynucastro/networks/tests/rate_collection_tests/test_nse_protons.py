import pynucastro as pyna


class TestNseProtons:

    def test_nse_protons(self, capsys, reaclib_library):

        nuc = ["p", "he4", "ti44", "v47", "cr48", "mn51",
               "fe52", "co55", "ni56"]

        lib = reaclib_library.linking_nuclei(nuc)
        rc = pyna.RateCollection(libraries=lib, verbose=True)

        rc.make_nse_protons(51)

        captured = capsys.readouterr()

        output = """modifying p_Mn51_to_Fe52_reaclib to use NSE protons
modifying Fe52_to_p_Mn51_reaclib to use NSE protons
modifying p_Co55_to_Ni56_reaclib to use NSE protons
modifying Ni56_to_p_Co55_reaclib to use NSE protons
modifying He4_Fe52_to_p_Co55_reaclib to use NSE protons
modifying p_Co55_to_He4_Fe52_reaclib to use NSE protons
"""

        assert captured.out == output
