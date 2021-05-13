#!/usr/bin/env python3

import numpy as np

from collections import namedtuple
from pynucastro import Composition, Nucleus
from pynucastro.rates import DummyNucleus
from pynucastro.nucdata import BindingTable
from mpi4py import MPI

def _wrap_conds(conds):
    
    try:
        conds[0]
        try:
            conds[0][0]
            return conds
        except:
            return [conds]
    except:
        raise ValueError('Conditions must be non-empty subscriptable object')

NetInfo = namedtuple("NetInfo", "y ydot z a ebind m")

def get_net_info(net, rho, T, comp):

    y_dict = comp.get_molar()
    ydots_dict = net.evaluate_ydots(rho, T, comp)

    bintable = BindingTable()
    bintable = {(nuc.n, nuc.z): nuc for nuc in bintable.nuclides}

    y = np.zeros(len(net.unique_nuclei), dtype=np.float64)
    ydot = np.zeros(len(net.unique_nuclei), dtype=np.float64)
    z = np.zeros(len(net.unique_nuclei), dtype=np.float64)
    a = np.zeros(len(net.unique_nuclei), dtype=np.float64)
    ebind = np.zeros(len(net.unique_nuclei), dtype=np.float64)
    m = np.zeros(len(net.unique_nuclei), dtype=np.float64)

    mass_neutron = 1.67492721184e-24
    mass_proton = 1.67262163783e-24
    c_light = 2.99792458e10

    for i, n in enumerate(net.unique_nuclei):

        y[i] = y_dict[n]
        ydot[i] = ydots_dict[n]
        z[i] = n.Z
        a[i] = n.A
        try:
            ebind[i] = bintable[n.N, n.Z].nucbind
        except KeyError:
            ebind[i] = 0.0
        m[i] = mass_proton * n.Z + mass_neutron * n.N - ebind[i] / c_light**2

    return NetInfo(y, ydot, z, a, ebind, m)

def enuc_dot(net_info):

    avo = 6.0221417930e23
    c_light = 2.99792458e10

    return -np.sum(net_info.ydot * net_info.m) * avo * c_light**2

def ye_dot(net_info):

    y, ydot, z, a, _, _ = net_info

    norm_fac = np.sum(y * a)
    return np.sum(ydot * z) / norm_fac - np.sum(y * z) / norm_fac**2 * np.sum(ydot * a)

def abar_dot(net_info):

    abar_inv = np.sum(net_info.y)
    print(1 / abar_inv)
    return -1 / abar_inv**2 * np.sum(net_info.ydot)

def map_comp1(comp, net):

    # Subtract out hydrogen, helium, carbon
    rem = 1. - comp.X[Nucleus('p')] - comp.X[Nucleus('he4')] \
             - comp.X[Nucleus('c12')]
    rem /= len(net.unique_nuclei) - 3

    omit = {Nucleus("p"), Nucleus("he4"), Nucleus("c12")}

    comp_new = Composition(net.unique_nuclei)

    for nuc in comp_new.X:
        if nuc in omit:
            comp_new.X[nuc] = comp.X[nuc]
        else:
            comp_new.X[nuc] = rem

    return comp_new
    
def map_comp2(comp, net):
    
    comp_new = Composition(net.unique_nuclei)

    for nuc in comp_new.X:
        comp_new.X[nuc] = comp.X[nuc]
        
    comp_new.normalize()
    return comp_new
    
def map_comp3(comp, net):
    
    comp_new = Composition(net.unique_nuclei)

    for nuc in comp_new.X:
        comp_new.X[nuc] = comp.X[nuc]
        
    return comp_new
    
def map_comp4(comp, net, dummy):
    
    comp_new = Composition(net.unique_nuclei)

    for nuc in comp_new.X:
        if nuc == dummy: continue
        comp_new.X[nuc] = comp.X[nuc]
        
    comp_new.X[dummy] = 1. - sum(comp_new.X.values())
    return comp_new

def rel_err(x, x0):

    return np.abs((x - x0) / x0)

def error_function(net_new, net_old, conds):

    conds = _wrap_conds(conds)
    err = np.zeros(3, dtype=np.float64)

    for rho, T, comp in conds:

        net_info1 = get_net_info(net_new, rho, T, comp)
        net_info2 = get_net_info(net_old, rho, T, comp)

        enuc_dot1 = enuc_dot(net_info1)
        enuc_dot2 = enuc_dot(net_info2)
        err[0] = max(err[0], rel_err(enuc_dot1, enuc_dot2))

        ye_dot1 = ye_dot(net_info1)
        ye_dot2 = ye_dot(net_info2)
        err[1] = max(err[1], rel_err(ye_dot1, ye_dot2))

        abar_dot1 = abar_dot(net_info1)
        abar_dot2 = abar_dot(net_info2)
        err[2] = max(err[2], rel_err(abar_dot1, abar_dot2))

    return err
    
def get_errfunc_enuc(net_old, conds):
    
    conds = _wrap_conds(conds)
    enucdot_list = []
    
    for rho, T, comp in conds:

        net_info_old = get_net_info(net_old, rho, T, comp)
        enucdot_list.append(enuc_dot(net_info_old))
    
    def erf(net_new):
        
        err = 0.0
        
        for cond, enucdot_old in zip(conds, enucdot_list):

            net_info_new = get_net_info(net_new, *cond)
            enucdot_new = enuc_dot(net_info_new)
            err = max(err, rel_err(enucdot_new, enucdot_old))
            
        return err
        
    return erf
    
def add_dummy_nucleus(red_net, conds):
    
    conds = _wrap_conds(conds)
    comps = (conds_i[2] for conds_i in conds)
    comp, = comps
        
    red_comp = map_comp3(comp, red_net)
    norm_fac = sum(red_comp.X.values())

    X_new = sum(comp.X.values()) - norm_fac
    Y_new = 1. / comp.eval_abar() - 1. / red_comp.eval_abar()

    A = X_new / Y_new
    Z = (comp.eval_ye() - norm_fac*red_comp.eval_ye()) / Y_new
    
    nuc = DummyNucleus(Z, A, 'dummy')
    red_net.add_nucleus(nuc)
    
    return nuc

if __name__ == "__main__":
    
    from pynucastro.reduction.load_network import load_network
    from pynucastro.reduction.generate_data import dataset
    from pynucastro.reduction import drgep, n_ary_search
    import time
    
    net = load_network(Nucleus('ni56'))
    data = list(dataset(net, n=6))
    
    # Perform DRGEP
    targets = map(Nucleus, ['p', 'ni56'])
    t0 = time.time()
    nuclei = drgep(net, data, targets, [1e-3, 1e-2], returnobj='nuclei')
    dt = time.time() - t0
    
    if MPI.COMM_WORLD.Get_rank() == 0:
        
        print()
        print(f"DRGEP reduction took {dt:.3f} s.")
        print("Number of species in full network: ", len(net.unique_nuclei))
        print("Number of species in DRGEP reduced network: ", len(nuclei))
        print()
        
    # Perform sensitivity analysis
    red_net = net.linking_nuclei(nuclei)
    errfunc = get_errfunc_enuc(net, data)
    t0 = time.time()
    red_net = n_ary_search(red_net, nuclei, errfunc)
    dt = time.time() - t0
    
    if MPI.COMM_WORLD.Get_rank() == 0:
        
        print(f"SA reduction took {dt:.3f} s.")
        print("Number of species in DRGEP + SA reduced network: ", len(red_net.unique_nuclei))
