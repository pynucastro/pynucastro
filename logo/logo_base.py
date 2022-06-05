import pynucastro as pyna

files = ["c12-pg-n13-ls09",
         "c13-pg-n14-nacr",
         "n13--c13-wc12",
         "n13-pg-o14-lg06",
         "n14-pg-o15-im05",
         "n15-pa-c12-nacr",
         "o14--n14-wc12",
         "o15--n15-wc12",
         "o14-ap-f17-Ha96c",
         "f17-pg-ne18-cb09",
         "ne18--f18-wc12",
         "f18-pa-o15-il10"]
rc = pyna.RateCollection(files)

comp = pyna.Composition(rc.get_nuclei())
comp.set_solar_like()

rc.plot(rho=2.e6, T=3.e7, comp=comp, outfile="cno_flow.png")
