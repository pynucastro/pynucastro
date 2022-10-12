#!/usr/bin/env python
"""
Uses AMETable class to extract binding energy from an AME Table.
"""
import argparse

from ame_table import AMETable

parser = argparse.ArgumentParser()
parser.add_argument('table', type=str, help='Name of AME Table file. E.g. "mass.mas12".')
parser.add_argument('-o', '--output', type=str, default='binding.txt',
                    help='Name of output file to store binding energies. Default is "binding.txt".')
args = parser.parse_args()

ame = AMETable(args.table)

fout = open(args.output, 'w')

fout.write('# AME Table file name: {}\n'.format(args.table))
fout.write('N                   Z                   Ebind(MeV)/Nucleon\n')

for nuc in ame.nuclides:
    ostr = '{:<20}'.format(nuc.n)
    ostr += '{:<20}'.format(nuc.z)
    ostr += '{:0.17e}\n'.format(nuc.nucbind)
    fout.write(ostr)

fout.close()
