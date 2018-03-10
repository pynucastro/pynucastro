# C-burning with A=23 URCA rate module generator

import pynucastro as pyna
from pynucastro.networks import StarKillerNetwork

library_file = "20180228default2"
mylibrary = pyna.rates.Library(library_file)

all_nuclei = ["he4", "c12", "o16", "n14", "f18", "ne21", "p", "n13", "ne20"]

subCh = mylibrary.linking_nuclei(all_nuclei, with_reverse=False)

net = StarKillerNetwork(libraries=subCh)
net.write_network()
