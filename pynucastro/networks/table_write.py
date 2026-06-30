import pynucastro as pyna
from pynucastro import Nucleus
import numpy as np
import nse_network as net
import nse_hdf5_network as net_hd

nucs = [Nucleus("p"), Nucleus("n"), Nucleus("he4"),
                Nucleus("fe52"), Nucleus("fe53"), Nucleus("fe54"),
                Nucleus("fe55"), Nucleus("fe56"),
                Nucleus("co54"), Nucleus("co55"), Nucleus("co56"),
                Nucleus("ni56"), Nucleus("ni57")]

tl = pyna.TabularLibrary().linking_nuclei(nucs)
rl = pyna.ReacLibLibrary().linking_nuclei(nucs)
all_lib = tl + rl

dupes = all_lib.find_duplicate_links()
rates_to_remove = []
for d in dupes:
    rates_to_remove += [r for r in d if isinstance(r, pyna.rates.ReacLibRate)]

for r in rates_to_remove:
    all_lib.remove_rate(r)

nse = net_hd.NSENetwork(libraries=[all_lib])

def get_reduced_comp(comp):
        nuc_list = [Nucleus("p"), Nucleus("n"), Nucleus("he4"),
                    Nucleus("fe52"), Nucleus("fe54"), Nucleus("ni56")]
        reduced_comp = comp.bin_as(nuc_list, exclude=[Nucleus("ni56")])
        X = []
        for n in reduced_comp.X:
            X.append((f"{n}", reduced_comp.X[n]))
        return X

Ts = np.logspace(9.6, 10.4, 5)
rhos = np.logspace(7, 10, 5)
yes = np.linspace(0.35, 0.5, 4)

nse.generate_table_hdf5(rho_values=rhos,
                   T_values=Ts,
                   Ye_values=yes,
                   comp_reduction_func=get_reduced_comp)

