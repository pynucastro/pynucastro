import re

with open('ffnrates.dat', 'r') as file:
    header = 0
    for l in file:

        s = l.strip()

        if s[0:3] == 'neg':
            l1 = l
            matches = re.findall(r"([a-zA-Z]{1,3}[0-9]{1,3})", l1)
            parent = matches[0]
            header += 1
            continue

        if s[0:3] == 'pos':
            l2 = l
            matches = re.findall(r"([a-zA-Z]{1,3}[0-9]{1,3})", l2)
            daughter = matches[0]
            header += 1
            continue

        if s[0:3] == 'end':
            continue

        if header == 2:
            with open(f'{parent}_{daughter}_table.dat', 'w') as output:
                output.write(l1)
                output.write(l2)
                output.write(l)
                header = 0
            continue

        with open(f'{parent}_{daughter}_table.dat', 'a') as output:
            output.write(l)
