# C-burning rate module generator

from pyreaclib.networks import BoxLibNetwork

files = ["c12-c12a-ne20-cf88",
         "c12-c12n-mg23-cf88",
         "c12-c12p-na23-cf88",
         "c12-ag-o16-nac2",
         "n--p-wc12"]

c_net = BoxLibNetwork(files, use_cse=True)
c_net.write_network()




