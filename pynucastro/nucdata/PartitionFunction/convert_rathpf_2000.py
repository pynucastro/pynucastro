#!/usr/bin/env python
"""
Reformat the partition function tables retrieved from https://download.nucastro.org/astro/fits/
into a standard format for pynucastro.

The tables part_frdm.asc and part_etfsiq.asc corresponds to:

Thomas Rauscher and Friedrich-Karl Thielemann, Atomic Data and Nuclear
Data Tables, 75:1–351, 2000
"""

import argparse

from pynucastro.nucdata import PeriodicTable

#Create a parser variable that demands for a table and creates an option to add a name to the formatted table
parser = argparse.ArgumentParser()
parser.add_argument('table', type=str, help='Name of the input table. "E.g. part_frdm.asc.txt"')
parser.add_argument('-o', '--output', type=str, default='frdm_low', help='Name of the formatted table E.g. ')

args = parser.parse_args()

finput = open(args.table, 'r')

temp = ['0.01', '0.15', '0.2', '0.3', '0.4', '0.5', '0.6', '0.7',
                '0.8', '0.9', '1.0', '1.5', '2.0', '2.5', '3.0', '3.5',
                '4.0', '4.5', '5.0', '6.0', '7.0', '8.0', '9.0', '10.0']

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

for _ in range(5):
    finput.readline()


for line in finput:

    nucleus_stats = line.strip().split()
    Z = int(nucleus_stats.pop(0))
    A = int(nucleus_stats.pop(0))
    J = float(nucleus_stats.pop(0))

    if Z == 0:
        assert A == 1
        abbrev = 'n'
    else:
        element = PeriodicTable.lookup_Z(Z)
        abbrev = '{}{}'.format(element.abbreviation, A)

    part1 = finput.readline().strip().strip('\n')
    part2 = finput.readline().strip().strip('\n')
    part3 = finput.readline().strip().strip('\n')

    fout.write('{}\n'.format(abbrev))
    fout.write('{}\n\n'.format(part1+' '+part2+' '+part3))

    finput.readline()

fout.close()
finput.close()
