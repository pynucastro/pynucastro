"""Extract the (A,Z, dm) tuples from `nubase_4.mas20.txt`, where:

:var A: is the atomic weight measured in atomic mass units.
:var Z: is the atomic number.
:var dm: is the mass difference A_{nuc}-A.

"""

import argparse

#os.path.dirname(os.path.realpath(__file__))

parser = argparse.ArgumentParser()
parser.add_argument('input', type=str, help='Name of the input table')
parser.add_argument('-o', '--output', type=str, default='mass_excess2020', help='Name of the formatted mass excess table')


args = parser.parse_args()

finput = open(args.input)

for _ in range(25):
    finput.readline()

fout = open(args.output+'.txt', 'w')
fout.write(f'# Mass difference evaluation table: {args.output} \n')
fout.write('# only ground states are tabulated \n')
fout.write('#\n')
fout.write('#\n')
fout.write('==A== {:18s} ==Z== {:10s} ======dm===== \n'.format(' ', ' '))

for line in finput:

    isomer_string = line[7]
    isomer = int(isomer_string)

    if isomer != 0:
        continue

    A_string = line[0:3].strip()
    Z_string = line[4:7].strip()
    dm_string = line[18:31].strip().strip('#')

    A = int(A_string)
    Z = int(Z_string)

    #dm is measured in keV, but we want MeV
    dm = float(dm_string)/1.0e3
    fout.write(f'{A:3d} {" ":20s} {Z:3d} {" ":10s} {dm:20.15} \n')

finput.close()
fout.close()
