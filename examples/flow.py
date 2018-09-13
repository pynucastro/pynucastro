from pynucastro.networks import RateCollection, Composition

files = ["c12-pg-n13-ls09",
         "c13-pg-n14-nacr",
         "n13--c13-wc12",
         "n13-pg-o14-lg06",
         "n14-pg-o15-im05",
         "n15-pa-c12-nacr",
         "o14--n14-wc12",
         "o15--n15-wc12"]
rc = RateCollection(files)

comp = Composition(rc.get_nuclei())

comp.set_all(0.005)
comp.set_nuc("p", 0.7)
comp.set_nuc("he4", 0.28)
comp.normalize()

rho = 2000
T = 2.e8

rr = rc.evaluate_rates(rho, T, comp)
for k, v in rr.items():
    print("{}: {}".format(k, v))


rc.plot(outfile="flow.png", rho=rho, T=T, comp=comp)
