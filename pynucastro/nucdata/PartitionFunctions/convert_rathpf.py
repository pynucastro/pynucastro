#!/usr/bin/env python
"""
Reformat the rathpf table retrieved from https://groups.nscl.msu.edu/jina/nucdatalib/
into a standard format for pynucastro.

This table corresponds to:

Thomas Rauscher and Friedrich-Karl Thielemann, Atomic Data and Nuclear
Data Tables, 75:1–351, 2000
"""
import argparse
from pynucastro.nucdata import PeriodicTable, Element

parser = argparse.ArgumentParser()
parser.add_argument('table', type=str,
                    help='Name of table file. E.g. "masslib_evaluation_rathpf.txt".')
parser.add_argument('-o', '--output', type=str, default='partition_functions_rathpf.txt',
                    help='Name of output file to store partition functions. Default is "partition_functions_rathpf.txt".')
args = parser.parse_args()

finput = open(args.table, 'r')

# Get rid of the first 13 lines
for _ in range(13):
    finput.readline()

# Get the list of temperature points in K
tline = finput.readline()
tsplit = tline.strip().split()
tsplit = tsplit[4:]
temperatures = ['{}e+9'.format(t) for t in tsplit]

# Throw away line 15
finput.readline()

fout = open(args.output, 'w')
fout.write('# Partition function evaluation name: rathpf\n')
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
        lsplit = line.split()
        Z = int(lsplit.pop(0))
        A = int(lsplit.pop(0))
        _ = lsplit.pop(0)
        pfun = lsplit
        if Z == 0:
            assert(A == 1)
            abbrev = 'n'
        else:
            enuc = PeriodicTable.lookup_Z(Z)
            abbrev = '{}{}'.format(enuc.abbreviation, A)
        fout.write('{}\n'.format(abbrev))
        fout.write('  '.join(pfun) + '\n\n')
fout.close()
finput.close()
