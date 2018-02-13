# CNO rate module generator
from pynucastro.networks import PythonNetwork

files = ["c12-pg-n13-ls09",
         "c13-pg-n14-nacr",
         "n13--c13-wc12",
         "n13-pg-o14-lg06",
         "n14-pg-o15-im05",
         "n15-pa-c12-nacr",
         "o14--n14-wc12",
         "o15--n15-wc12"]

cno_net = PythonNetwork(files)
cno_net.write_network("cno_rhs.py")
