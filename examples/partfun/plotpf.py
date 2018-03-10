import numpy as np
import matplotlib.pyplot as plt
from pynucastro.nucdata import PartitionFunction, PartitionFunctionCollection
from pynucastro.rates import Nucleus

import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-nuc', '--nucleus', type=str, help="Name of nucleus for which to plot the partition function.")
parser.add_argument('-etfsiq', action="store_true", help="If supplied, use the ETFSIQ high temperature partition functions from Rauscher 2003.")
parser.add_argument('-frdm', action="store_true", help="If supplied, use the FRDM high temperature partition functions from Rauscher 2003.")
args = parser.parse_args()

def plot_pf(pf):
    t = pf.temperature
    f = pf.partition_function
    g = pf.interpolant

    xt = np.linspace(t[0], t[-1], num=1000)

    plt.loglog(t, f, 'r.', label=pf.nucleus)
    plt.loglog(xt, np.power(10.0, g(xt/1.0e9)), 'b-', label='spline k = 3')
    plt.xlabel('Log T (K)')
    plt.ylabel('Log F')
    plt.savefig('partition_function_{}.png'.format(pf.nucleus), dpi=600)
    plt.clf()

if __name__ == "__main__":
    pfc = PartitionFunctionCollection()
    if args.etfsiq and args.frdm:
        print('Cannot use both ETFSIQ and FRDM high temperature partition functions.')
        exit()
    high_temp_fun = False
    if args.etfsiq:
        high_temp_fun = "rauscher2003_ETFSIQ"
    if args.frdm:
        high_temp_fun = "rauscher2003_FRDM"
    pfc.use_high_temperature_functions(high_temp_fun)

    if args.nucleus:
        nucleus = Nucleus(args.nucleus)
        pf = pfc.get_partition_function(nucleus, high_temperature_partition_functions=high_temp_fun)
        plot_pf(pf)
    else:
        for pf in pfc:
            plot_pf(pf)
