import pynucastro as pyna

lib = pyna.ReacLibLibrary().linking_nuclei(["he4", "c12", "o16", "ne20", "na20",
                                            "n", "p", "na23", "ne23", "mg23", "mg24"])
rates_to_derive = lib.backward().get_rates()
for r in rates_to_derive:
    fr = lib.get_rate_by_nuclei(r.products, r.reactants)
    if fr:
        lib.remove_rate(r)
        d = pyna.DerivedRate(fr, use_pf=False)
        lib.add_rate(d)
net = pyna.FortranNetwork(libraries=[lib])
fig = net.plot(rotated=True)
fig.savefig("test.png")
net.write_network()
