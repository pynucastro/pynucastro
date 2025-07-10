"""Extract the (A,Z, tau) tuples from `nubase_4.mas20.txt`, where:

:var A: is the atomic weight measured in atomic mass units.
:var Z: is the atomic number.
:var tau: is the halflife in s (or `"stable"` if it is stable)

"""

import argparse

parser = argparse.ArgumentParser()
parser.add_argument('input', type=str, help='Name of the input table')
parser.add_argument('-o', '--output', type=str, default='halflife2020', help='Name of the formatted mass excess table')

seconds_per_year = 31556926

units = {"ms": 1.e-3,
         "us": 1.e-6,
         "ns": 1.e-9,
         "ps": 1.e-12,
         "fs": 1.e-15,
         "as": 1.e-18,
         "zs": 1.e-21,
         "ys": 1.e-24,
         "s": 1,
         "m": 60,
         "h": 3600,
         "d": 24*3600,
         "y": seconds_per_year,
         "ky": 1.e3 * seconds_per_year,
         "My": 1.e6 * seconds_per_year,
         "Gy": 1.e9 * seconds_per_year,
         "Ty": 1.e12 * seconds_per_year,
         "Py": 1.e15 * seconds_per_year,
         "Ey": 1.e18 * seconds_per_year,
         "Zy": 1.e21 * seconds_per_year,
         "Yy": 1.e24 * seconds_per_year}

args = parser.parse_args()

finput = open(args.input)

for _ in range(25):
    finput.readline()

fout = open(args.output+'.txt', 'w')
fout.write(f'# Halflife evaluation table: {args.output} \n')
fout.write('#\n')
fout.write('#\n')
fout.write('#\n')
fout.write('==A== {:18s} ==Z== {:10s} ======tau===== \n'.format(' ', ' '))

for line in finput:

    # skip isomers
    isomer_string = line[7]
    isomer = int(isomer_string)

    if isomer != 0:
        continue

    tau = line[69:78].strip().strip("#")

    # sometimes it is a lower limit, like ">5 zs" -- we'll just
    # strip off the ">"
    tau = tau.replace(">", "")
    tau = tau.replace("<", "")
    tau = tau.replace("~", "")

    tau_unit = line[78:80].strip()

    print(tau, tau_unit)
    if tau == "stbl":
        tau = "stable"
    elif tau_unit in units:
        tau = float(tau) * units[tau_unit]
    else:
        continue

    A_string = line[0:3].strip()
    Z_string = line[4:7].strip()

    A = int(A_string)
    Z = int(Z_string)

    if tau == "stable":
        fout.write(f'{A:3d} {" ":20s} {Z:3d} {" ":10s} {tau:25.16} \n')
    else:
        fout.write(f'{A:3d} {" ":20s} {Z:3d} {" ":10s} {tau:25.16g} \n')

finput.close()
fout.close()
