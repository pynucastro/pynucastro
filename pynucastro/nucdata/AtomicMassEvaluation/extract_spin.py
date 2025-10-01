"""Extract the spin data for each ground state (gs) nuclei,
characterized by the pair (A, Z), from nubase_3.mas20.txt, published
in:

Kondev, F. G., Wang, M., Huang, W. J., Naimi, S., & Audi, G.
Chinese Physics C, 45(3), 030001. (2021) doi:10.1088/1674-1137/abddae

located in Table I.

"""

import argparse


def num_states(spin_str_element):
    """Evaluate the spin number string, formatted as s=a/b and
    returns the number of states 2*s + 1.

    In the table we have three type of strings:
    1. spin numbers integers formatted with 1 or 2 characters, e.g s=1, and s=10.
    2. spin numbers formatted with 3 characters. e.g. s=3/2.
    3. spin numbers formatted with 4 characters. e.g. s=11/2

    Parameters
    ----------
    :var spin_str_element: This string class element contains the information
                      in [88:102], about spin, parity, and
                      isospin charge.

    Return
    ------
    :var states: this integer variable contains the number of states associated
            to the `spin_str_element` string

    """

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
parser.add_argument('-o', '--output', type=str, default='spins2020', help='Pynucastro Formatted Table')

args = parser.parse_args()

finput = open(args.table, 'r')

# First, we need to get rid of the first 25 lines.
for _ in range(25):
    finput.readline()

fout = open(args.output+'.txt', 'w')

# Write the content of the new formatted types of A, Z, Spin, Number of States
# in a new file.
fout.write('# Ground state spin evaluation table: {}\n'.format(args.output))
fout.write('# if the ==Spin== column contain more than two values, the gs state is uncertain.\n')
fout.write('#\n')
fout.write('#==A=='+' '*22+'==Z=='+' '*18+'==Spin=='+' '*18+'==Number=of=States=='+' '*12+'==Experimental=='+' '*12+'==Theoretical==\n')

# For each consecutive line, we extract:
for line in finput:
    A_string = line[0:3]                               # The string that contains A
    Z_string = line[4:7]                               # The string that contains Z
    i_str = line[7]                                    # The string that defines gs
    spin_str_list = line[88:102].strip().split()       # A list of strings ["spin+parity", "Isospin"]

    # We convert the first three string variables to integers
    A = int(A_string)
    Z = int(Z_string)
    i = int(i_str)

    # We extract the "spin+parity" string from spin_str_list
    if spin_str_list:
        spin_str = spin_str_list.pop(0)
    else:
        continue

    # We remove the characters from the string "+-,*()#".
    # The meaning of each (set of) character(s) is(are) the following:
    #
    # "*"   : The spin/parity measurement is provided by
    #         strong experimental arguments.
    # "+-"  : The parity associated to each nucleus.
    # "()"  : The spin/parity measurement is provided by
    #         weak experimental arguments.
    # "#"   : The spin/parity measurement is provided by
    #         theoretical arguments.

    experimental_str = ' '
    theoretical_str = ' '

    for c in spin_str:
        if c == '*':
            experimental_str = 's'
            theoretical_str = 'w'
        elif c == '#':
            experimental_str = 'w'
            theoretical_str = 's'
        elif c != '(' and c != ')' and c != '*' and c != '#':
            experimental_str = 's'
            theoretical_str = 's'

    if spin_str[0] == '(':
        experimental_str = 'w'
        theoretical_str = 'w'

    special_chars = "+-*()#,"

    for c in special_chars:
        spin_str = spin_str.replace(c, ' ')

    spin_str = spin_str.strip().split()

    if i == 0:
        # We eliminate lines that may contain more than one spin value in the gs
        if len(spin_str) > 1:
            continue

        output_str = '{0:5}'.format(A)

        output_str += ' '*21
        output_str += '{0:5}'.format(Z)

        output_str += ' '*23
        spin = spin_str.pop(0)
        output_str += '{0:<4}'.format(spin)
        states = str(num_states(spin))

        output_str += ' '*27
        output_str += '{0:<3}'.format(states)

        output_str += ' '*27
        output_str += '{0:<3}'.format(experimental_str)

        output_str += ' '*27
        output_str += '{0:<3}'.format(theoretical_str)

        output_str += '\n'
        fout.write(output_str)
    else:
        continue

fout.close()
finput.close()
