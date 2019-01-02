# C-burning with A=23 URCA rate module generator

from pynucastro.networks import StarKillerNetwork

files = ["c12-c12a-ne20-cf88",
         "c12-c12n-mg23-cf88",
         "c12-c12p-na23-cf88",
         "c12-ag-o16-nac2",
         "na23--ne23-toki",
         "ne23--na23-toki",
         "n--p-wc12"]

urca_net = StarKillerNetwork(files)
urca_net.write_network()
