'''
This program plots the abundances in a file produced by intnet.

Headers should include 'Time' and nuclides of the form 'Y_{nuc}',
where {nuc} is, e.g. he4, c12, ...

Example:
$ python plot_abundances.py net_history.dat

Donald Willcox
'''
from cycler import cycler
import matplotlib.pyplot as plt
import numpy as np
import re
import argparse

parser = argparse.ArgumentParser()
# Take the net_history input file as a command line argument
parser.add_argument('nethistoryfile', type=str)
args = parser.parse_args()

f = open(args.nethistoryfile,'r')
h = f.readline().split()
d = {}
for hi in h:
        d[hi] = []

for l in f:
    ls = l.split()
    for hi, lsi in zip(h, ls):
        d[hi].append(float(lsi))

for hi in h:
    d[hi] = np.array(d[hi])

for hi in h:
    if hi != 'Time':
        fig = plt.figure()
        ax=fig.add_subplot(111)
        ax.plot(d['Time'], d[hi])
        ax.set_yscale('log')
        ax.set_xscale('log')
        plt.xlabel('T (s)')
        his = hi.split('_')
        hip = his[0] + '_{' + his[1] + '}'
        plt.ylabel('$\\mathrm{'+hip+'}$')
        plt.savefig(hi + '.pdf')

# Setup plot setting cycling
plt.rc('lines', linewidth=2)
plt.rc('axes', prop_cycle=(cycler('color', ['r','g','b','y']) * cycler('linestyle', ['-', '--',':','-.'])))

# Plot mass fractions X = Y*A
fig = plt.figure()
ax=fig.add_subplot(111)
legend_handles = []
for hi in h:
        if 'Y_' in hi:
                his = hi.split('_')
                shortname = his[1]
                m = re.search('[a-zA-Z]*([0-9]*)',shortname)
                if m.group(1):
                        A = int(m.group(1))
                else:
                        A = 1
                hip = 'X' + '_{' + shortname + '}'
                ax.plot(d['Time'], d[hi]*A, label='$'+hip+'$')
ax.set_xlim([100.0,22589.6685427])
ax.set_ylim([1.0e-3,1.0])
ax.set_xscale('log')
ax.set_yscale('log')
plt.xlabel('T (s)')
framepos = ax.get_position()
ax.set_position([framepos.x0, framepos.y0, framepos.width*0.8, framepos.height])
plt.legend(bbox_to_anchor=(1.4,1))
plt.ylabel('X')
plt.savefig('net_history.pdf')
