# triple-alpha rate module generator
import pyreaclib

files = ["c12-gaa-he4-fy05",
         "he4-aag-c12-fy05"]

pyreaclib.make_network(files, "triple-alpha_rhs.py")




