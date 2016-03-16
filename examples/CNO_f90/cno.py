# CNO rate module generator

import reaclib

files = ["c12-pg-n13-ls09",
         "n13--c13-wc12",
         "c13-pg-n14-nacr",
         "n14-pg-o15-im05",
         "o15--n15-wc12",
         "n15-pa-c12-nacr",
         "n13-pg-o14-lg06",
         "o14--n14-wc12"]


rc = reaclib.RateCollection(files)

rc.make_network_f90()




