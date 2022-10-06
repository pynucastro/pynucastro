#!/usr/bin/env python
"""
Reformat Table 2 and Table 3 of Rauscher 2003
into a standard format for pynucastro.

These tables corresponds to:

Thomas Rauscher, Astrophysical Journal Supplement Series, 147:403-408, 2003

These table txt files are published together with the paper above
and are read unmodified by this script.
"""

import argparse

from pynucastro.nucdata import PeriodicTable

#Create a parser variable that demands for a table and creates an option to add a name to the formatted table
parser = argparse.ArgumentParser()
parser.add_argument('table', type=str, help='Name of the input table. "E.g. datafile2.txt"')
parser.add_argument('-o', '--output', type=str, default='frdm_high', help='Name of the formatted table e.g fdrm_high.txt')

args = parser.parse_args()

finput = open(args.table, 'r')

temp = ['12.0', '14.0', '16.0', '18.0', '20.0', '22.0', '24.0', '26.0',
        '28.0', '30.0', '35.0', '40.0', '45.0', '50.0', '55.0', '60.0',
        '65.0', '70.0', '75.0', '80.0', '85.0', '90.0', '95.0', '100.0',
        '105.0', '110.0', '115.0', '120.0', '125.0', '130.0', '135.0', '140.0',
        '145.0', '150.0', '155.0', '160.0', '165.0', '170.0', '175.0', '180.0',
        '190.0', '200.0', '210.0', '220.0', '230.0', '240.0', '250.0', '275.0']

temperatures = ['{}E+9'.format(t) for t in temp]

fout = open(args.output+'.txt', 'w')
fout.write('# Partition function evaluation name: {}\n'.format(args.output))
fout.write('# Each entry is of the form:\n')
fout.write('#\n')
fout.write('# [nucleus name, e.g. ni56]\n')
fout.write('# [List of partition function values]\n')
fout.write('#\n')
fout.write('# Partition functions all evaluated at these temperatures (K):\n')
fout.write('  '.join(temperatures) + '\n\n')

for _ in range(63):
    finput.readline()

count = 0
for line in finput:
    count += 1
    nucleus_stats = line.strip().split()[1:]
    Z = int(nucleus_stats.pop(0))
    A = int(nucleus_stats.pop(0))
    J = float(nucleus_stats.pop(0))

    if Z == 0:
        assert A == 1
        abbrev = 'n'
    else:
        element = PeriodicTable.lookup_Z(Z)
        abbrev = '{}{}'.format(element.abbreviation, A)

    fout.write('{}\n'.format(abbrev))
    fout.write(' '.join(nucleus_stats) + '\n\n')

fout.close()
finput.close()
