# Test module for generating a network containing
# one rate from each Reaclib chapter

from pynucastro.networks import StarKillerNetwork

files = ["b17-nnn-c14-wc12",
         "he3-he3pp-he4-nacr",
         "he4-aag-c12-fy05",
         "he4-npahe3-li7-mafo",
         "he4-pphe3-he3-nacr",
         "he6-gnn-he4-cf88",
         "li7-tnna-he4-mafo",
         "n--p-wc12",
         "p-ng-d-an06",
         "t-gn-d-nk06",
         "t-pn-he3-de04"]

net = StarKillerNetwork(files)
net.write_network()
