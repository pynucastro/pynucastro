'''
Figure out what the rho-T grid looks like for data tables
from Toki, et al. 2015.

Donald Willcox
'''

import numpy as np
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('infile',type=str, help='The input file to process.')
args = parser.parse_args()

try: ifile = open(args.infile, 'r')
except: raise

# Eat header
for n in xrange(0,6):
    ifile.readline()

dens = []
temp = []

for l in ifile:
    if l.strip() != '':
        ls = l.split()
        dens.append(ls[0])
        temp.append(ls[1])
ifile.close()

dens = list(set(dens))
temp = list(set(temp))

dens_f = [float(s) for s in dens]
temp_f = [float(s) for s in temp]

dens = np.array(dens_f)
temp = np.array(temp_f)

print(args.infile)
print('')
print('dens, #  : ' + str(len(dens)))
print('dens, min: ' + str(np.amin(dens)))
print('dens, max: ' + str(np.amax(dens)))
print('')
print('temp, #  : ' + str(len(temp)))
print('temp, min: ' + str(np.amin(temp)))
print('temp, max: ' + str(np.amax(temp)))
