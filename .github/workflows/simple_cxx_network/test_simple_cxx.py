import pynucastro as pyna

lib = pyna.ReacLibLibrary().linking_nuclei(["he4", "c12", "o16", "ne20", "na20"])
net = pyna.SimpleCxxNetwork(libraries=[lib])
net.write_network()
