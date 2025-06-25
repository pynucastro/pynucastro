import pynucastro as pyna

lib = pyna.ReacLibLibrary().linking_nuclei(["he4", "c12", "o16", "ne20", "na20",
                                            "n", "p", "na23", "ne23", "mg23", "mg24"])
net = pyna.FortranNetwork(libraries=[lib])
fig = net.plot(rotated=True)
fig.savefig("test.png")
net.write_network()
