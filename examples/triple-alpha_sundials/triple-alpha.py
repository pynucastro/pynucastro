# triple-alpha rate module generator for a fortran network

from pynucastro.networks import SundialsNetwork

files = ["c12-gaa-he4-fy05",
         "he4-aag-c12-fy05"]

triple_alpha_net = SundialsNetwork(files)
triple_alpha_net.write_network()




