"""
Test table lookup algorithm for f90 implementation

Donald Willcox
"""

import numpy as np
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('infile',type=str, help='The input file to process.')
args = parser.parse_args()

try: ifile = open(args.infile, 'r')
except: raise

num_rhoy = 152
num_temp = 39
num_header = 7

# Eat header
for n in xrange(num_header):
    ifile.readline()

d = {}
d['logrhoy'] = []
d['logtemp'] = []
d['data'] = np.empty((num_temp, num_rhoy, 6), dtype=float)

jtab_mu = 0
jtab_dq = 1
jtab_vs = 2
jtab_rate = 3
jtab_nuloss = 4
jtab_gamma = 5

for j in xrange(num_rhoy):
    for i in xrange(num_temp):
        l = ifile.readline()
        # print j
        # print i
        # print l
        ls = l.split()
        d['logrhoy'].append(ls[0])
        d['logtemp'].append(ls[1])
        d['data'][i,j,jtab_mu] = float(ls[2])
        d['data'][i,j,jtab_dq] = float(ls[3])
        d['data'][i,j,jtab_vs] = float(ls[4])
        d['data'][i,j,jtab_rate] = float(ls[5])
        d['data'][i,j,jtab_nuloss] = float(ls[6])
        d['data'][i,j,jtab_gamma] = float(ls[7])
    if (j != num_rhoy-1):
        ifile.readline()
ifile.close()

d['logrhoy'] = np.array(sorted([float(s) for s in list(set(d['logrhoy']))]))
d['logtemp'] = np.array(sorted([float(s) for s in list(set(d['logtemp']))]))

def vector_index_lu(vector, fvar):
    print('IN VECTOR_INDEX_LU')
    print(fvar)
    print(vector)
    index = 0   # return value
    n = len(vector)-1
    if ( fvar < vector[0] ):
       index = 0
       return index
    elif ( fvar > vector[n] ):
       index = n
       return index
    else:
       nup = n
       ndn = 0
       for i in xrange(n+1):
          j = ndn + (nup - ndn)/2
          if ( fvar < vector[j] ):
             nup = j
          else:
             ndn = j
          if ( ((nup - ndn) == 1) ):
             index = ndn
             print(index)
             return index

def bilinear_lookup(rhoy, temp, ivar):
    fvar = 0.0  # return value
    logrhoy = np.log10(rhoy)
    logtemp = np.log10(temp)

    # Get box-corner points for interpolation
    # This deals with out-of-range inputs via linear extrapolation
    irhoy_lo = vector_index_lu(d['logrhoy'], logrhoy)
    itemp_lo = vector_index_lu(d['logtemp'], logtemp)
    irhoy_hi = irhoy_lo + 1
    itemp_hi = itemp_lo + 1

    print('irhoy_lo:')
    print(irhoy_lo)
    print('irhoy_hi:')
    print(irhoy_hi)
    print('itemp_lo:')
    print(itemp_lo)
    print('itemp_hi:')
    print(itemp_hi)

    # Bilinear interpolation within the box (log=log10)
    # T ^   B .      . C
    #   |
    # g |  AB   ABCD   CD
    # o |     .      .
    # l |   A          D
    #   |___________________> log rho*Ye
    temp_lo = d['logtemp'][ itemp_lo ]
    temp_hi = d['logtemp'][ itemp_lo+1 ]
    rhoy_lo = d['logrhoy'][ irhoy_lo ]
    rhoy_hi = d['logrhoy'][ irhoy_lo+1 ]

    print('rhoy_lo:')
    print(rhoy_lo)
    print('rhoy_hi:')
    print(rhoy_hi)
    print('temp_lo:')
    print(temp_lo)
    print('temp_hi:')
    print(temp_hi)

    fa = d['data'][ itemp_lo, irhoy_lo, ivar ]
    fb = d['data'][ itemp_hi, irhoy_lo, ivar ]
    fc = d['data'][ itemp_hi, irhoy_hi, ivar ]
    fd = d['data'][ itemp_lo, irhoy_hi, ivar ]
    fab = ( fa * ( temp_hi - logtemp ) + fb * ( logtemp - temp_lo ) )/( temp_hi - temp_lo )
    fcd = ( fd * ( temp_hi - logtemp ) + fc * ( logtemp - temp_lo ) )/( temp_hi - temp_lo )
    fvar = ( fab * ( rhoy_hi - logrhoy ) + fcd * ( logrhoy - rhoy_lo ) )/(rhoy_hi - rhoy_lo)
    return fvar

if __name__ == "__main__":
    while(True):
        rhoy = float(raw_input('Enter a log(rho*ye) value:'))
        temp = float(raw_input('Enter a log(temp) value:'))
        ivar = int(raw_input('Enter the value index to lookup (rate index is 3):'))
        print('Bilinear interpolation: ')
        print(bilinear_lookup(10.0**rhoy, 10.0**temp, ivar))
