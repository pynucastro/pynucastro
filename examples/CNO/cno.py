# CNO rate module generator

import reaclib

files = ["c12-pg-n13-ls09", 
         "c13-pg-n14-nacr",
         "n13--c13-wc12",
         "n13-pg-o14-lg06",
         "n14-pg-o15-im05",
         "n15-pa-c12-nacr",
         "o14--n14-wc12",
         "o15--n15-wc12"]

rf = open("cno_rates.py", "w")

rf.write("import numpy as np\n")

for f in files:
    r = reaclib.Rate(f)
    rf.write(r.function_string())

rf.close()

             



