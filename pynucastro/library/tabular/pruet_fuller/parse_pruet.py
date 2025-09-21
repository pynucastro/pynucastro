import re

import numpy as np

import pynucastro as pyna

with open("pruet_fowler_datafile2.txt") as df:

    # skip the header

    for _ in range(31):
        df.readline()

    rate_str = re.compile(r"Z=(\d+),N=(\d+)-->Z=(\d+),N=(\d+)")

    # now read and break into groups by what reaction it is

    current_reaction = ""

    ecf = None
    betaf = None

    while line := df.readline():
        T = np.log10(float(line[0:6]) * 1.e9)
        rhoYe = line[7:11]
        uf = line[16:22]

        # electron capture, beta+-decay
        r1 = line[23:48]
        rate = rate_str.match(r1)

        Zp_ec = int(rate.group(1))
        Np_ec = int(rate.group(2))
        ec_parent = pyna.Nucleus.from_Z_A(Zp_ec, Zp_ec + Np_ec)

        Zd_ec = int(rate.group(3))
        Nd_ec = int(rate.group(4))
        ec_daughter = pyna.Nucleus.from_Z_A(Zd_ec, Zd_ec + Nd_ec)

        r1em = float(line[49:57])
        r1cap = float(line[60:68])
        log_ec_rate = np.log10(10.0**r1em + 10.0**r1cap)
        log_ec_nu = line[71:79]

        # beta decay, positron capture
        r2 = line[80:105]
        rate = rate_str.match(r2)

        Zp_beta = int(rate.group(1))
        Np_beta = int(rate.group(2))
        beta_parent = pyna.Nucleus.from_Z_A(Zp_beta, Zp_beta + Np_beta)

        Zd_beta = int(rate.group(3))
        Nd_beta = int(rate.group(4))
        beta_daughter = pyna.Nucleus.from_Z_A(Zd_beta, Zd_beta + Nd_beta)

        r2em = float(line[106:114])
        r2cap = float(line[117:125])
        log_beta_rate = np.log10(10.0**r2em + 10.0**r2cap)
        log_beta_nu = line[128:136]

        if r1 != current_reaction:
            # we are dealing with a new set of rates, so open
            # new files

            current_reaction = r1

            # close existing files
            if ecf:
                ecf.close()
            if betaf:
                betaf.close()

            # open the new files

        # output

        print(T, rhoYe, uf, r1, r1em, r1cap, log_ec_nu, r2, r2em, r2cap, log_beta_nu)
        print(ec_parent, ec_daughter, beta_parent, beta_daughter)
        break
