import argparse
import os
from pynucastro.nucdata import PeriodicTable, Element

def num_states(spin_str_element):

    if len(spin_str_element) == 1 or len(spin_str_element) == 2:
        states = 2*int(spin_str_element) + 1
        return states

    elif len(spin_str_element) == 3:
        num = int(spin_str_element[0])
        den = int(spin_str_element[2])
        states = 2*num//den + 1
        return states

    elif len(spin_str_element) == 4:
        num = int(spin_str_element[0:2])
        den = int(spin_str_element[3])
        states = 2*num//den + 1
        return states

    else:
        return None

parser = argparse.ArgumentParser()
parser.add_argument('table', type=str, help='Name of the input spin stable')
parser.add_argument('-o', '--output', type=str, default = 'nubase2020', help='Pynucastro Formatted Table')

args = parser.parse_args()

finput = open(args.table, 'r')
for _ in range(25):
    finput.readline()

fout = open(args.output+'.txt', 'w')

fout.write('# Ground state spin evaluation table: {}\n'.format(args.output))
fout.write('# if the ==Spin== column contain more than two values, the gs state is uncertain.\n')
fout.write('#\n')
fout.write('#==A=='+' '*22+'==Z=='+' '*18+'==Spin=='+' '*18+'==Number=of=States==\n')

for line in finput:
    A_string = line[0:3] 
    Z_string = line[4:7]
    i_str = line[7] 
    spin_str_list = line[88:102].strip().split()

    A = int(A_string)
    Z = int(Z_string)
    i = int(i_str)

    if spin_str_list:
        spin_str = spin_str_list.pop(0)
    else:
        continue

    special_chars = "+-*()#,"

    for c in special_chars:
        spin_str  = spin_str.replace(c, ' ')

    spin_str = spin_str.strip().split()

    if i == 0:
        output_str = '{0:5}'.format(A)

        output_str += ' '*21
        output_str += '{0:5}'.format(Z)

        if spin_str:
            output_str += ' '*23
            spin_str_1 = spin_str.pop(0)
            output_str += '{0:<4}'.format(spin_str_1) 
            state1 = str(num_states(spin_str_1))
            
            if spin_str:
                output_str += ' '
                spin_str_2 = spin_str.pop(0)
                output_str += '{0:<5}'.format(spin_str_2)
                state2 = str(num_states(spin_str_2))
            else:
                state2 = None
        else:
            continue

        if state2:
            output_str += ' '*21
        else: 
            output_str +=' '*27

        output_str += '{0:<3}'.format(state1)

        if state2:
            output_str += ' ' 
            output_str += '{0:<3}'.format(state2)

        output_str += '\n'
        fout.write(output_str)
    else:
        continue

fout.close()
finput.close() 


