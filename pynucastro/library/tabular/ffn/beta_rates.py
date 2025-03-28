import os
import re

import numpy as np
import pandas as pd
import scipy.constants as scp

pd.set_option('display.max_colwidth', 100)

# Input Files:
#=============

curr_dir = os.getcwd()
tables = [f.name for f in os.scandir(curr_dir) if re.search("_table.dat", f.name)]

for loc in tables:

    # Extracting Q and nuclides:
    #===========================
    with open(loc, 'r') as file:

        for l in file:
            s = l.strip()

            if s[0:3] == 'neg':
                l1 = l
                matches = re.findall(r"([a-zA-Z]{1,3}[0-9]{1,3})", l1)
                parent = matches[0]
                parent_g = re.match(r"([A-Za-z]*)([0-9]*)", parent)
                parent_str = f"{parent_g[2]}{parent_g[1]}"
                qmatch = l1.split('Q=').pop()
                Qcap = float(qmatch)
                continue

            if s[0:3] == 'pos':
                l2 = l
                matches = re.findall(r"([a-zA-Z]{1,3}[0-9]{1,3})", l2)
                daughter = matches[0]
                daughter_g = re.match(r"([A-Za-z]*)([0-9]*)", daughter)
                daughter_str = f"{daughter_g[2]}{daughter_g[1]}"
                qmatch = l2.split('Q=').pop()
                Qbeta = float(qmatch)
                continue

    # Defining Output:
    #=================

    with open(f"{parent_str}-{daughter_str}_electroncapture.dat", 'w') as outfile:
        outfile.write(f"!{parent_str} -> {daughter_str}, e- capture\n")
        outfile.write(f"!Q={Qcap} MeV\n")
        outfile.write("!\n")

    with open(f"{daughter_str}-{parent_str}_betadecay.dat", 'w') as outfile:
        outfile.write(f"!{daughter_str} -> {parent_str}, beta-decay\n")
        outfile.write(f"!Q={Qbeta} MeV\n")
        outfile.write("!\n")

    MeV_to_eV = 1.0e6
    eV_to_J, _, _ = scp.physical_constants['electron volt-joule relationship']
    J_to_erg = 1.0e7

    MeV_to_erg = MeV_to_eV * eV_to_J * J_to_erg

    #Loading the DataFrame:
    #======================
    #Fields: 't9', 'lrho', 'uf', 'lbeta+', 'leps-', 'lnu', 'lbeta-', 'leps+', 'lnubar'
    table = pd.read_table(loc, sep=r'\s+', header=0, comment='#', skiprows=3)
    table = table.dropna()

    #Common Variables:
    #=================
    df = {}
    df['!Log(rhoY)'] = table['lrho']
    df['Log(temp)'] = np.log10(table['t9']) + 9.0
    df['mu'] = table['uf'] * MeV_to_erg
    df['dQ'] = 0.00
    df['Vs'] = 0.00

    #Electron-Capture
    #================
    df_p = pd.DataFrame.from_dict(df)
    df_p['!Log(rhoY)'] = df_p['!Log(rhoY)'].apply(lambda x: f"{float(x):.6f}")
    df_p['Log(temp)'] = df_p['Log(temp)'].apply(lambda x: f"{float(x):.6f}")
    df_p['mu'] = df_p['mu'].apply(lambda x: f"{float(x):.6e}")
    df_p['dQ'] = df_p['dQ'].apply(lambda x: f"{float(x):.2f}")
    df_p['Vs'] = df_p['Vs'].apply(lambda x: f"{float(x):.2f}")

    df_p['Log(e-cap-rate)'] = np.log10(10**table['lbeta+'] + 10**table['leps-'])
    df_p['Log(e-cap-rate)'] = df_p['Log(e-cap-rate)'].apply(lambda x: f"{float(x):.6f}")
    df_p['Log(nu-energy-loss)'] = table['lnubar'] + np.log10(MeV_to_erg)
    df_p['Log(nu-energy-loss)'] = df_p['Log(nu-energy-loss)'].apply(lambda x: f"{float(x):.6f}")
    df_p['Log(gamma-energy)'] = -100.00
    df_p['Log(gamma-energy)'] = df_p['Log(gamma-energy)'].apply(lambda x: f"{float(x):.2f}")

    df_p.loc[-1] = ['!Log(g/cm^3)', 'Log(K)', 'erg', 'erg', 'erg', 'Log(1/s)', 'Log(erg/s)', 'Log(erg/s)']  # units
    df_p.index = df_p.index + 1  # shifting index
    df_p = df_p.sort_index()  # sorting by index
    df_p.loc[0] = df_p.loc[0].apply(lambda x: f"{x:<15s}")

    #Beta Decay
    #==========
    df_m = pd.DataFrame.from_dict(df)
    df_m['!Log(rhoY)'] = df_m['!Log(rhoY)'].apply(lambda x: f"{float(x):.6f}")
    df_m['Log(temp)'] = df_m['Log(temp)'].apply(lambda x: f"{float(x):.6f}")
    df_m['mu'] = df_m['mu'].apply(lambda x: f"{float(x):.6e}")
    df_m['dQ'] = df_m['dQ'].apply(lambda x: f"{float(x):.2f}")
    df_m['Vs'] = df_m['Vs'].apply(lambda x: f"{float(x):.2f}")

    df_m['Log(beta-decay-rate)'] = np.log10(10**table['lbeta-'] + 10**table['leps+'])
    df_m['Log(beta-decay-rate)'] = df_m['Log(beta-decay-rate)'].apply(lambda x: f"{float(x):.6f}")
    df_m['Log(nu-energy-loss)'] = table['lnubar'] + np.log10(MeV_to_erg)
    df_m['Log(nu-energy-loss)'] = df_m['Log(nu-energy-loss)'].apply(lambda x: f"{float(x):.6f}")
    df_m['Log(gamma-energy)'] = -100.00
    df_m['Log(gamma-energy)'] = df_m['Log(gamma-energy)'].apply(lambda x: f"{float(x):.2f}")

    df_m.loc[-1] = ['!Log(g/cm^3)', 'Log(K)', 'erg', 'erg', 'erg', 'Log(1/s)', 'Log(erg/s)', 'Log(erg/s)']  # units
    df_m.index = df_p.index + 1  # shifting index
    df_m = df_p.sort_index()  # sorting by index

    df_m.loc[0] = df_m.loc[0].apply(lambda x: f"{x:<15s}")

    # Appending tables and writing output toki files:
    # ==============================================

    with open(f"{parent_str}-{daughter_str}_electroncapture.dat", 'a') as outfile:
        df_p.to_string(outfile, justify='left', index=False)

    with open(f"{daughter_str}-{parent_str}_betadecay.dat", 'a') as outfile:
        df_m.to_string(outfile, justify='left', index=False)

    with open(f"{parent}--{daughter}-toki", 'w') as outfile:
        outfile.write("t\n")
        outfile.write(f"       {parent} {daughter}\n")
        outfile.write(f"{parent_str}-{daughter_str}_electroncapture.dat\n")
        outfile.write("5\n")
        outfile.write("11\n")
        outfile.write("13\n")

    with open(f"{daughter}--{parent}-toki", 'w') as outfile:
        outfile.write("t\n")
        outfile.write(f"       {daughter} {parent}\n")
        outfile.write(f"{daughter_str}-{parent_str}_beta_decay.dat\n")
        outfile.write("5\n")
        outfile.write("11\n")
        outfile.write("13\n")
