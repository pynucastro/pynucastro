# triple-alpha rate module generator

from pynucastro.networks import PythonNetwork

files = ["c12-gaa-he4-fy05",
         "he4-aag-c12-fy05"]

triple_alpha_net = PythonNetwork(files)
triple_alpha_net.write_network("triple-alpha_rhs.py")
