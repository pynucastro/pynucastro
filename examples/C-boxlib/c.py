# C-burning rate module generator

import reaclib

files = ["c12-c12a-ne20-cf88",
         "c12-c12n-mg23-cf88",
         "c12-c12p-na23-cf88",
         "c12-ag-o16-nac2",
         "n--p-wc12"]

rc = reaclib.RateCollection(files)

rc.make_network('boxlib')




