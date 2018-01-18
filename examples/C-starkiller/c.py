# C-burning rate module generator

from pyreaclib.networks import StarKillerNetwork

files = ["c12-c12a-ne20-cf88",
         "c12-c12n-mg23-cf88",
         "c12-c12p-na23-cf88",
         "c12-ag-o16-nac2",
         "n--p-wc12"]

c_net = StarKillerNetwork(files)
c_net.write_network(use_cse=True)
