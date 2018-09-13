"""
This program eats 7 lines of header and reformats the lines in a
specified input file so numbers can be compared with
the same formatting between the input file and a reference (not supplied).

Donald Willcox
"""

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

lines = []
for l in ifile:
    if l.strip() != '':
        ls = l.split()
        lines.append(ls)
ifile.close()

lines_n = []

for l in lines:
    lt = []
    for s in l:
        lt.append(float(s))
    lines_n.append(lt)

try: of = open(args.infile+'.compare','w')
except: raise
for lt in lines_n:
    for lti in lt:
        of.write('{: >25.14e}'.format(lti).replace('e','E'))
    of.write('\n')

of.close()

try: of = open(args.infile+'.compare.min','w')
except: raise
for lt in lines_n:
    for lti in lt[2:]:
        of.write('{: >25.14e}'.format(lti).replace('e','E'))
    of.write('\n')

of.close()
