#!/usr/bin/env python3

import argparse

parser = argparse.ArgumentParser()
parser.add_argument("infile", type=str)
parser.add_argument("-p", action="store_true", help="print stuff for checking which rate is which")
args = parser.parse_args()

table_file = None
header_lines = None


def nuc_A_first(nuc_A_last):
    snuc = ""
    for c in nuc_A_last:
        if c.isdigit():
            snuc += c
    for c in nuc_A_last:
        if not c.isdigit():
            snuc += c
    return snuc


try:
    # get the nuclei from the name of the rate metadata file
    nucs_file = args.infile[:-5].split("--")
    f = open(args.infile, "r")
    chapter = f.readline().strip()
    # assert the chapter is "t" for a tabulated rate
    assert chapter == "t"
    # get the nuclei inside the rate metadata file
    nucs = f.readline().strip().split()
    # assert the nuclei in the rate metadata match the nuclei in the metadata filename
    assert nucs_file[0] == nucs[0]
    assert nucs_file[1] == nucs[1]
    # get the name of table file and number of header lines
    table_file = f.readline().strip()
    # get the nuclei from the table file and assert
    # those are the same nuclei involved in the rate
    table_file_nucs = table_file.split("_")[0].split("-")
    assert table_file_nucs[0].lower() == nuc_A_first(nucs[0]).lower()
    assert table_file_nucs[1].lower() == nuc_A_first(nucs[1]).lower()
    header_lines = int(f.readline().strip())
    # assert table sizes are correct
    n = int(f.readline().strip())
    assert n == 152
    n = int(f.readline().strip())
    assert n == 39
    f.close()

    if args.p:
        print(args.infile)
        print(table_file)

    # checks to make sure valid data begins after the
    # number of header lines the rate metadata file says
    f = open(table_file, "r")
    for i, line in enumerate(f):
        if i > header_lines:
            break
        if i < header_lines:
            ls = line.strip()
            assert ls == "" or ls.startswith("!")
        if i == header_lines:
            lss = line.strip().split()
            rhoy = lss[0]
            temp = lss[1]
            assert rhoy == "1.000000e+07"
            assert temp == "1.000000e+07"
        if args.p:
            print("{}: {}".format(i, line.strip()))
    f.close()
except:  # noqa
    print("FAILED")
else:
    print("OK")
