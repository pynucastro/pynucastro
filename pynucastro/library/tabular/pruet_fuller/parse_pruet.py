# create a class that holds all the data for a line
# create a dict keyed by (ec_parent, ec_daughter)
# apprend to a list in each key
# loop over keys, find unique Ts

import re

import numpy as np

import pynucastro as pyna
from pynucastro.constants import constants

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

        # the column numbers used in the slicing come from the header
        # in the original data table

        T = np.log10(float(line[0:6]) * 1.e9)
        rhoYe = line[7:11]
        uf = float(line[16:22]) * constants.MeV2erg

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
        # sum the emission and capture and convert nu loss from MeV/s to erg/s
        # a value of -100 means it was not computed originally
        if r1em == -100.0 and r1cap == -100.0:
            log_ec_rate = -100.0
        else:
            log_ec_rate = np.log10(10.0**r1em + 10.0**r1cap)
        nu = float(line[71:79])
        if nu == -100.0:
            log_ec_nu = -100.0
        else:
            log_ec_nu = np.log10(10.0**nu * constants.MeV2erg)

        # beta decay, positron capture
        # the rate here might not exist
        beta_valid = True

        r2 = line[80:105]
        rate = rate_str.match(r2)

        Zp_beta = int(rate.group(1))
        Np_beta = int(rate.group(2))
        beta_parent = pyna.Nucleus.from_Z_A(Zp_beta, Zp_beta + Np_beta)

        Zd_beta = int(rate.group(3))
        Nd_beta = int(rate.group(4))
        beta_daughter = pyna.Nucleus.from_Z_A(Zd_beta, Zd_beta + Nd_beta)

        try:
            r2em = float(line[106:114])
            r2cap = float(line[117:125])
            # sum the emission and capture and convert nu loss from MeV/s to erg/s
            # a value of -100 means it was not computed originally
            if r2em == -100.0 and r2cap == -100.0:
                log_beta_rate = -100.0
            else:
                log_beta_rate = np.log10(10.0**r2em + 10.0**r2cap)
            nu = float(line[128:136])
            if nu == -100.0:
                log_beta_nu = -100.0
            else:
                log_beta_nu = np.log10(10.0**nu * constants.MeV2erg)
        except ValueError:
            beta_valid = False

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
            ecf = open(f"pruetfowler-{ec_parent.A}{ec_parent.el}-{ec_daughter.A}{ec_daughter.el}_electroncapture.dat", "w")
            if beta_valid:
                betaf = open(f"pruetfowler-{beta_parent.A}{beta_parent.el}-{beta_daughter.A}{beta_daughter.el}_betadecay.dat", "w")
            else:
                betaf = None

            # write the headers
            ecf.write(f"!{ec_parent.A}{ec_parent.el} -> {ec_daughter.A}{ec_daughter.el}, e- capture\n")
            ecf.write(f"!Q={ec_daughter.mass - ec_parent.mass:5.3f} MeV\n")
            ecf.write("!\n")
            ecf.write(f"!{'Log(rhoY)':20} {'Log(temp)':20} {'mu':20} {'dQ':20} {'Vs':20} {'Log(e-cap-rate)':20} {'Log(nu-energy-loss)':20} {'Log(gamma-energy)':20}\n")
            ecf.write(f"!{'Log(g/cm^3)':20} {'Log(K)':20} {'erg':20} {'erg':20} {'erg':20} {'Log(1/s)':20} {'Log(erg/s)':20} {'Log(erg/s)':20}\n")

            if betaf:
                betaf.write(f"!{beta_parent.A}{beta_parent.el} -> {beta_daughter.A}{beta_daughter.el}, beta-decay\n")
                betaf.write(f"!Q={beta_daughter.mass - beta_parent.mass:5.3f} MeV\n")
                betaf.write("!\n")
                betaf.write(f"!{'Log(rhoY)':20} {'Log(temp)':20} {'mu':20} {'dQ':20} {'Vs':20} {'Log(beta-decay-rate)':20} {'Log(nu-energy-loss)':20} {'Log(gamma-energy)':20}\n")
                betaf.write(f"!{'Log(g/cm^3)':20} {'Log(K)':20} {'erg':20} {'erg':20} {'erg':20} {'Log(1/s)':20} {'Log(erg/s)':20} {'Log(erg/s)':20}\n")

        # output

        ecf.write(f" {rhoYe:20} {T:10.7f} {uf:20.14g} {0.0:20} {0.0:20} {log_ec_rate:20} {log_ec_nu:20} {-100.0:20}\n")
        if betaf:
            betaf.write(f" {rhoYe:20} {T:10.7f} {uf:20.14g} {0.0:20} {0.0:20} {log_beta_rate:20} {log_beta_nu:20} {-100.0:20}\n")
