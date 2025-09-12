import matplotlib.pyplot as plt

plt.rcParams['mathtext.fontset'] = 'custom'
plt.rcParams['mathtext.rm'] = 'Fira Sans'
plt.rcParams['mathtext.it'] = 'Fira Sans:italic'
plt.rcParams['mathtext.bf'] = 'Fira Sans:bold'


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

fig = rc.plot(rho=2.e6, T=3.e7, comp=comp,
              node_size=1300, node_font_size=18, node_color="#444444")

# ax[0] is the plot, ax[1] is the colorbar
ax = fig.get_axes()
ax[0].set_xlim(5.5, 8.5)
ax[0].set_ylim(5.5, 8.2)

print(ax)

ax[0].set_axis_off()
ax[1].remove()

fig.text(0.42, 0.1, "pynucastro", horizontalalignment="center",
         color="#444444",
         fontsize="40", fontstyle="italic", fontfamily="Fira Sans")

ax[0].set_aspect("equal", adjustable="datalim")

fig.set_size_inches(4.5, 4.0)


ax[0].margins(0.0)
fig.tight_layout()

fig.savefig("logo.pdf", bbox_inches="tight", pad_inches=0)
fig.savefig("logo.svg", bbox_inches="tight", pad_inches=0)
fig.savefig("logo.png", bbox_inches="tight", pad_inches=0)
