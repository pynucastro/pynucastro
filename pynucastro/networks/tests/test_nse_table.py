import sys
from pathlib import Path

import numpy as np
import pytest

import pynucastro as pyna
from pynucastro import Nucleus


@pytest.mark.skipif(sys.platform == "darwin")
class TestNSETable:
    @pytest.fixture(scope="class")
    def nse_net(self, reaclib_library, tabular_library):

        nucs = [Nucleus("p"), Nucleus("n"), Nucleus("he4"),
                Nucleus("fe52"), Nucleus("fe53"), Nucleus("fe54"),
                Nucleus("fe55"), Nucleus("fe56"),
                Nucleus("co54"), Nucleus("co55"), Nucleus("co56"),
                Nucleus("ni56"), Nucleus("ni57")]

        tlib = tabular_library.linking_nuclei(nucs, print_warning=False)
        rlib = reaclib_library.linking_nuclei(nucs)

        all_lib = rlib + tlib

        dupes = all_lib.find_duplicate_links()

        rates_to_remove = []
        for d in dupes:
            rates_to_remove += [r for r in d if isinstance(r, pyna.rates.ReacLibRate)]

        for r in rates_to_remove:
            all_lib.remove_rate(r)

        nse = pyna.NSENetwork(libraries=[all_lib])

        return nse

    @staticmethod
    def get_reduced_comp(comp):
        nuc_list = [Nucleus("p"), Nucleus("n"), Nucleus("he4"),
                    Nucleus("fe52"), Nucleus("fe54"), Nucleus("ni56")]
        reduced_comp = comp.bin_as(nuc_list, exclude=[Nucleus("ni56")])
        X = []
        for n in reduced_comp.X:
            X.append((f"{n}", reduced_comp.X[n]))
        return X

    def test_generate_table(self, nse_net):

        Ts = np.logspace(9.6, 10.4, 3)
        rhos = np.logspace(7, 10, 4)
        yes = np.linspace(0.43, 0.5, 3)

        nse_net.generate_table(rho_values=rhos,
                               T_values=Ts,
                               Ye_values=yes,
                               comp_reduction_func=self.get_reduced_comp)

        # this creates a file called `nse.tbl`, which we want to compare
        # to the stored benchmark

        base_path = Path(__file__).parent.relative_to(Path.cwd())
        ref_path = base_path/"_nse_table"

        with (
            open("nse.tbl") as new_table,
            open(f"{ref_path}/nse_scipy_1.14.tbl") as ref_table_1,
            open(f"{ref_path}/nse_scipy_1.15.tbl") as ref_table_2
        ):

            new_lines = new_table.readlines()
            ref_lines_1 = ref_table_1.readlines()
            ref_lines_2 = ref_table_2.readlines()

            for new, ref_1, ref_2 in zip(new_lines, ref_lines_1, ref_lines_2):
                if new.startswith("#"):
                    continue
                # pylint: disable-next=consider-using-in
                assert new == ref_1 or new == ref_2
