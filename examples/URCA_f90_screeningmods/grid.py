"""
Integrate a grid of dens and temp using the defaults in net.par.default

Donald Willcox
"""
import numpy as np
import subprocess
import re

num_dens = 100
num_temp = 100
dens = np.logspace(8.0, 11.0, num=num_dens, endpoint=True)
temp = np.logspace(7.0, 9.65, num=num_temp, endpoint=True)
enuc_ave = np.zeros((100,100))

def fmt_to_dp(f):
    return '{:1.14e}'.format(float(f)).replace('e', 'd')

def set_dp_namelist(vard):
    # vard is a dictionary keyed by var names with var values
    try:
        fi = open('net.par.default', 'r')
    except:
        raise
    try:
        fo = open('net.par', 'w')
    except:
        raise
    red = {}
    for k in vard.keys():
        red[k] = re.compile('[\s]*'+k+'[\s]*=.*')
    for line in fi:
        foundre = False
        for k in vard.keys():
            rem = red[k].match(line)
            if rem:
                oline = '  '+k+' = '+fmt_to_dp(vard[k])+'\n'
                foundre = True
        if not foundre:
            oline = line
        fo.write(oline)
    fi.close()
    fo.close()

renuc = re.compile('Average enuc_dot:[\s]+([0-9\.]+)[dDEe](\+?\-?[0-9]*)')
for i, rho in enumerate(dens):
    for j, t in enumerate(temp):
        set_dp_namelist({'cv_pars%dens':rho,
                         'cv_pars%temp':t})
        comm = subprocess.Popen(['./integrator.gfortran.exe'], shell=False,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE)
        comm_out, comm_err = comm.communicate()
        if comm_err:
            print('ERROR , idens={} , dens={} , itemp={} , temp={}'.format(i,rho,j,t))
            print(comm_err)
            exit()
        comm_out_utf = comm_out.decode('utf-8').split('\n')
        for line in comm_out_utf:
            rem = renuc.match(line)
            if rem:
                enucave = float(rem.group(1)) * (10.0**int(rem.group(2)))
                print('Found Average enuc_dot: {}'.format(enucave))
                enuc_ave[i][j] =  enucave

# Save grid data
try:
    fgrid = open('grid.dat', 'w')
except:
    raise
fgrid.write("!num dens: {}\n".format(num_dens))
fgrid.write("!num temp: {}\n".format(num_temp))
fgrid.write("! dens                 temp                 enuc_dot\n")
for i, rho in enumerate(dens):
    for j, t in enumerate(temp):
        fgrid.write(' '.join([fmt_to_dp(rho),
                              fmt_to_dp(t),
                              fmt_to_dp(enuc_ave[i][j])])+'\n')
fgrid.close()





