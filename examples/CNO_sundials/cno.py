# CNO rate module generator
from pynucastro.networks import SundialsNetwork

files = ["c12-pg-n13-ls09",
         "n13--c13-wc12",
         "c13-pg-n14-nacr",
         "n14-pg-o15-im05",
         "o15--n15-wc12",
         "n15-pa-c12-nacr",
         "n13-pg-o14-lg06",
         "o14--n14-wc12"]

cno_net = SundialsNetwork(files)
cno_net.write_network()




