"""
This program plots the emission and capture tables produced
by the program output_table.f90 akin to the plots in Toki, et al 2013.

Donald Willcox
"""

import numpy as np
import argparse
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument('--emission_infile',type=str, help='The emission input file to plot.')
parser.add_argument('--capture_infile',type=str, help='The capture input file to plot.')
args = parser.parse_args()

try: emt_file = open(args.emission_infile, 'r')
except: raise

try: cap_file = open(args.capture_infile, 'r')
except: raise

def read_rate_file(infile):
    ## Process header
    d = []

    # Get the number of temperature and rho*ye points
    n_temp = int(infile.readline().strip())
    n_rhoy = int(infile.readline().strip())

    # Eat header
    infile.readline()

    # Get data
    for i in xrange(n_temp):
        dt = {}
        dt['log_lrhoy'] = []
        dt['log_ltemp'] = 0.0
        dt['rate'] = []
        for j in xrange(n_rhoy):
            l = infile.readline().strip().split()
            dt['log_lrhoy'].append(float(l[0]))
            dt['log_ltemp'] = float(l[1])
            dt['rate'].append(float(l[2]))
        dt['log_lrhoy'] = np.array(dt['log_lrhoy'])
        dt['rate'] = np.array(dt['rate'])
        d.append(dt)
    return d

d_emt = read_rate_file(emt_file)
d_cap = read_rate_file(cap_file)
emt_file.close()
cap_file.close()

## Plot emission and capture rates
fig = plt.figure()
ax=fig.add_subplot(111)
for dt in d_emt:
    ax.plot(dt['log_lrhoy'], dt['rate'], '-.')
for dt in d_cap:
    ax.plot(dt['log_lrhoy'], dt['rate'], '-')
plt.xlabel('$\\mathrm{Log_{10}~\\rho Y_e~[g~cm^{-3}]}$')
plt.ylabel('$\\mathrm{Log_{10}~\\lambda~[s^{-1}]}$')
ax.set_ylim([-25,0])
plt.tight_layout()
plt.title('A=23')
plt.savefig('output_table_rates.pdf')
