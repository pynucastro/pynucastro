import re

with open('odarates.dat', 'r') as file:
    header = 0
    for l in file:

        s = l.strip()

        if s[0:3] == 'pos':
            l1 = l
            matches = re.findall(r"([a-zA-Z]{1,3}\s[0-9]{1,3})", l1)
            daughter = matches[0]
            header += 1
            continue

        if s[0:2] == 'T9':
            if header == 0:
                continue
            header += 1
            l2 = l
            continue

        if header == 2:
            with open(f'{daughter}_table.dat', 'w') as output:
                output.write(l1)
                output.write(l2)
                output.write(l)
                header = 0
            continue

        with open(f'{daughter}_table.dat', 'a') as output:
            output.write(l)
