# CNO rate module generator
import pynucastro as pyna
from pynucastro.networks import PythonNetwork

library_file = "20180201ReaclibV2.22"
mylibrary = pyna.rates.Library(library_file)

all_nuclei = ["p", "he4", "c12", "n13", "c13", "n14", "c13", "o14", "o15", "n15"]

cno_lib = mylibrary.linking_nuclei(all_nuclei)

cno_net = PythonNetwork(libraries=cno_lib)
cno_net.write_network("network_module.py")
