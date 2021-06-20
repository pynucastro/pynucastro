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

fout.write(f'# AME Table file name: {args.table}\n')
fout.write('N                   Z                   Ebind(MeV)/Nucleon\n')

for nuc in ame.nuclides:
    ostr =  f'{nuc.n:<20}'
    ostr += f'{nuc.z:<20}'
    ostr += f'{nuc.nucbind:0.17e}\n'
    fout.write(ostr)

fout.close()
