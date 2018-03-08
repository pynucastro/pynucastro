#!/usr/bin/env python
"""
Reformat Table 2 and Table 3 of Rauscher 2003
into a standard format for pynucastro.

This table corresponds to:

Thomas Rauscher, Astrophysical Journal Supplement Series, 147:403-408, 2003

These table txt files are published together with the paper above 
and are read unmodified by this script.
"""
import argparse
from pynucastro.nucdata import PeriodicTable, Element

parser = argparse.ArgumentParser()
parser.add_argument('table', type=str,
                    help='Name of table file. E.g. "rauscher2003.txt".')
parser.add_argument('-o', '--output', type=str, default='partition_functions_rauscher2003.txt',
                    help='Name of output file to store partition functions. Default is "partition_functions_rauscher2003.txt".')
args = parser.parse_args()

finput = open(args.table, 'r')

# Get rid of the first 4 lines
for _ in range(4):
    finput.readline()

# Get the 5th line
f5 = finput.readline()

# Determine whether this is the FRDM or ETFSIQ dataset
if 'FRDM' in f5:
    name = 'rauscher2003_FRDM'
else:
    name = 'rauscher2003_ETFSIQ'

# Get rid of the next 9 lines
for _ in range(9):
    finput.readline()

# Get the list of temperature points in K
temperatures = []
for line in finput:
    if '-----' in line:
        break
    else:
        lsplit = line.strip().split('T = ')
        try:
            _ = float(lsplit[-1])
        except:
            print('Error -- could not read temperature {}'.format(lsplit[-1]))
        else:
            temperatures.append(lsplit[-1])

fout = open(args.output, 'w')
fout.write('# Partition function evaluation name: {}\n'.format(name))
fout.write('# Each entry is of the form:\n')
fout.write('#\n')
fout.write('# [nucleus name, e.g. ni56]\n')
fout.write('# [List of partition function values]\n')
fout.write('#\n')
fout.write('# Partition functions all evaluated at these temperatures (K):\n')
fout.write('  '.join(temperatures) + '\n\n')

# Read nucleus partition functions
for line in finput:
    if line.strip():
        lsplit = line.strip().split()
        abbrev = lsplit.pop(0)
        pfun   = lsplit[3:]
        fout.write('{}\n'.format(abbrev))
        fout.write('  '.join(pfun) + '\n\n')
fout.close()
finput.close()
