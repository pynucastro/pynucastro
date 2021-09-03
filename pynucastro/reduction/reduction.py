#!/usr/bin/env python3

import numpy as np

from collections import namedtuple
from pynucastro import Composition, Nucleus
from pynucastro.rates import DummyNucleus
from pynucastro.nucdata import BindingTable
import path_flux_analysis as pfa
from parallel_pfa import ppfa
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

def get_net_info(net, comp, rho, T):

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

def get_net_info_arr(net, rho, T, comp, s_p, s_c):

    y_dict = comp.get_molar()
    ydot = net.evaluate_ydots_arr(rho, T, comp, s_p, s_c)

    bintable = BindingTable()
    bintable = {(nuc.n, nuc.z): nuc for nuc in bintable.nuclides}

    y = np.zeros(len(net.unique_nuclei), dtype=np.float64)
    z = np.zeros(len(net.unique_nuclei), dtype=np.float64)
    a = np.zeros(len(net.unique_nuclei), dtype=np.float64)
    ebind = np.zeros(len(net.unique_nuclei), dtype=np.float64)
    m = np.zeros(len(net.unique_nuclei), dtype=np.float64)

    mass_neutron = 1.67492721184e-24
    mass_proton = 1.67262163783e-24
    c_light = 2.99792458e10

    for i, n in enumerate(net.unique_nuclei):

        y[i] = y_dict[n]
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
        
    enucdot_list = []
    
    for comp, rho, T in conds:
        net_info_old = get_net_info(net_old, comp, rho, T)
        enucdot_list.append(enuc_dot(net_info_old))
    
    def erf(net_new):
        
        err = 0.0
        
        for cond, enucdot_old in zip(conds, enucdot_list):

            net_info_new = get_net_info(net_new, *cond)
            enucdot_new = enuc_dot(net_info_new)
            err = max(err, rel_err(enucdot_new, enucdot_old))
            
        return err
        
    return erf
    
def get_errfunc_enuc_mpi(net_old, conds):
    
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    N_proc = comm.Get_size()

    #only designing this to work for 2^n processes, 2^m conditions for m >= 3
    n = np.array([len(conds[0]), len(conds[1]), len(conds[2])])
    if(n[2] >= N_proc):
        comp_i = (n[2]//N_proc)*(rank)
        comp_f = (n[2]//N_proc)*(rank+1)
        rho_i = 0
        rho_f = n[0]
    else:
        comp_split = (N_proc // n[2])
        comp_i = rank // (comp_split)
        comp_f = comp_i + 1
        rho_i = (n[0]//comp_split) * (rank % comp_split)
        rho_f = (n[0]//comp_split) * ((rank+1) % comp_split)
        if rank % comp_split == comp_split-1:
            rho_f = n[0]

    rho_L = conds[0][rho_i:rho_f]
    T_L = conds[1]
    comp_L = conds[2][comp_i:comp_f]
    
    n_loc = np.array([len(comp_L), len(rho_L), len(T_L)])

    n_map, r_map = pfa.get_maps(net_old)
    s_p, s_c, s_a = pfa.get_stoich_matrices(net_old, r_map)

    enucdot_list = []
    for comp in comp_L:
        net_old.update_yfac_arr(composition=comp, s_c=s_c)
        for rho in rho_L:
            net_old.update_prefac_arr(rho=rho, composition=comp)
            for T in T_L:
                net_info_old = get_net_info_arr(net_old, rho, T, comp, s_p, s_c)
                enucdot_list.append(enuc_dot(net_info_old))
    
    # enucdot_list_g = np.zeros(N_proc*len(enucdot_list))
    # comm.Allgather([np.array(enucdot_list), MPI.DOUBLE], [enucdot_list_g, MPI.DOUBLE])

    def erf(net_new):
        
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        N_proc = comm.Get_size()
        
        err = 0.0

        n_map, r_map = pfa.get_maps(net_new)
        s_p, s_c, s_a = pfa.get_stoich_matrices(net_new, r_map)
        net_new.update_coef_arr()

        for i, comp in enumerate(comp_L):
            comp = map_comp3(comp, net_new)
            net_new.update_yfac_arr(composition=comp, s_c=s_c)
            for j, rho in enumerate(rho_L):
                net_new.update_prefac_arr(rho=rho, composition=comp)
                for k, T in enumerate(T_L):
                    idx = k + j*n_loc[2] + i*n_loc[1]*n_loc[2]
                    enucdot_old = enucdot_list[idx]
                    net_info_new = get_net_info_arr(net_new, rho, T, comp, s_p, s_c)
                    enucdot_new = enuc_dot(net_info_new)
                    err = max(err, rel_err(enucdot_new, enucdot_old))
        
        err_arr_loc = np.array([err])
        err_arr = np.zeros_like(err_arr_loc)
        comm.Allreduce([err_arr_loc, MPI.DOUBLE], [err_arr, MPI.DOUBLE], op=MPI.MAX)
            
        return err_arr[0]
        
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
    from pynucastro.reduction import drgep, binary_search_trim, sens_analysis
    import time
    import sys

    if len(sys.argv) == 5:
        endpoint = Nucleus(sys.argv[1])
        n = (sys.argv[2], sys.argv[3], sys.argv[4])
    elif len(sys.argv) == 2:
        endpoint = Nucleus(sys.argv[1])
        n = 16
    elif len(sys.argv) == 1:
        endpoint = Nucleus('te108')
        n = 16
    else:
        print("Usage: ./load_network.py <endpoint>")
    
    net = load_network(endpoint)
    data = list(dataset(net, n=n))
    
    # Perform DRGEP / PFA
    targets = map(Nucleus, ['p', 'ni56'])
    t0 = time.time()
    nuclei = drgep(net, data, targets, [1e-3, 1e-2], returnobj='nuclei', use_mpi=True)
    dt = time.time() - t0
    
    if MPI.COMM_WORLD.Get_rank() == 0:
        
        print()
        print(f"DRGEP reduction took {dt:.3f} s for {MPI.COMM_WORLD.Get_size()} processes.")
        print("Number of species in full network: ", len(net.unique_nuclei))
        print("Number of species in DRGEP reduced network: ", len(nuclei))
        print()
    
    # Perform sensitivity analysis
    red_net = net.linking_nuclei(nuclei)
    errfunc = get_errfunc_enuc(net, data)
    MPI.COMM_WORLD.Barrier()
    t0 = time.time()
    # red_net = binary_search_trim(red_net, nuclei, errfunc)
    red_net = sens_analysis(red_net, errfunc, 0.05, True)
    dt = time.time() - t0
    
    if MPI.COMM_WORLD.Get_rank() == 0:
    
        print(f"SA reduction took {dt:.3f} s.")
        print("Number of species in DRGEP + SA reduced network: ", len(red_net.unique_nuclei))
        print("Error: ", errfunc(red_net))
