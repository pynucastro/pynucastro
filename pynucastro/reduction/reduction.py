#!/usr/bin/env python3

import numpy as np

from collections import namedtuple
from pynucastro import Composition, Nucleus
from pynucastro.nucdata import DummyNucleus, BindingTable
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
            ebind[i] = bintable.get_nuclide(n.N, n.Z).nucbind
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
    
# def get_errfunc_enuc_mpi(net_old, conds):
# 
#     comm = MPI.COMM_WORLD
#     rank = comm.Get_rank()
#     N_proc = comm.Get_size()
# 
#     #only designing this to work for 2^n processes, 2^m conditions for m >= 3
#     n = np.array([len(conds[0]), len(conds[1]), len(conds[2])])
#     if(n[2] >= N_proc):
#         comp_i = (n[2]//N_proc)*(rank)
#         comp_f = (n[2]//N_proc)*(rank+1)
#         rho_i = 0
#         rho_f = n[0]
#     else:
#         comp_split = (N_proc // n[2])
#         comp_i = rank // (comp_split)
#         comp_f = comp_i + 1
#         rho_i = (n[0]//comp_split) * (rank % comp_split)
#         rho_f = (n[0]//comp_split) * ((rank+1) % comp_split)
#         if rank % comp_split == comp_split-1:
#             rho_f = n[0]
# 
#     rho_L = conds[0][rho_i:rho_f]
#     T_L = conds[1]
#     comp_L = conds[2][comp_i:comp_f]
# 
#     n_loc = np.array([len(comp_L), len(rho_L), len(T_L)])
# 
#     n_map, r_map = pfa.get_maps(net_old)
#     s_p, s_c, s_a = pfa.get_stoich_matrices(net_old, r_map)
# 
#     enucdot_list = []
#     for comp in comp_L:
#         net_old.update_yfac_arr(composition=comp, s_c=s_c)
#         for rho in rho_L:
#             net_old.update_prefac_arr(rho=rho, composition=comp)
#             for T in T_L:
#                 net_info_old = get_net_info_arr(net_old, rho, T, comp, s_p, s_c)
#                 enucdot_list.append(enuc_dot(net_info_old))
# 
#     # enucdot_list_g = np.zeros(N_proc*len(enucdot_list))
#     # comm.Allgather([np.array(enucdot_list), MPI.DOUBLE], [enucdot_list_g, MPI.DOUBLE])
# 
#     def erf(net_new):
# 
#         comm = MPI.COMM_WORLD
#         rank = comm.Get_rank()
#         N_proc = comm.Get_size()
# 
#         err = 0.0
# 
#         n_map, r_map = pfa.get_maps(net_new)
#         s_p, s_c, s_a = pfa.get_stoich_matrices(net_new, r_map)
#         net_new.update_coef_arr()
# 
#         for i, comp in enumerate(comp_L):
#             comp = map_comp3(comp, net_new)
#             net_new.update_yfac_arr(composition=comp, s_c=s_c)
#             for j, rho in enumerate(rho_L):
#                 net_new.update_prefac_arr(rho=rho, composition=comp)
#                 for k, T in enumerate(T_L):
#                     idx = k + j*n_loc[2] + i*n_loc[1]*n_loc[2]
#                     enucdot_old = enucdot_list[idx]
#                     net_info_new = get_net_info_arr(net_new, rho, T, comp, s_p, s_c)
#                     enucdot_new = enuc_dot(net_info_new)
#                     err = max(err, rel_err(enucdot_new, enucdot_old))
# 
#         err_arr_loc = np.array([err])
#         err_arr = np.zeros_like(err_arr_loc)
#         comm.Allreduce([err_arr_loc, MPI.DOUBLE], [err_arr, MPI.DOUBLE], op=MPI.MAX)
# 
#         return err_arr[0]
# 
#     return erf
    
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
    red_net.add_inert_nucleus(nuc)
    
    return nuc

if __name__ == "__main__":
    
    import time
    import sys
    import argparse
    from load_network import load_network
    from generate_data import dataset
    from pynucastro.reduction import drgep, binary_search_trim, sens_analysis, pfa
    
    #-----------------------------------
    # Setup parser and process arguments
    #-----------------------------------
    
    description = "Example/test script for running a selected reduction algorithm."
    algorithm_help = "The algorithm to use. Currently supports 'drgep' and 'pfa'."
    endpoint_help = "The nucleus to use as an endpoint in the unreduced network."
    datadim_help = """The dimensions of the dataset to generate. Ordering is the number of densities,
            then temperatures, then metallicities. Default is [4, 4, 4]."""
    drho_help = "Number of density points to include in the dataset."
    dtemp_help = "Number of temperature points to include in the dataset."
    dmetal_help = "Number of metallicity points to include in the dataset."
    brho_help = "Range of densities to include in the dataset (in g/cm^3)."
    btemp_help = "Range of temperatures to include in the dataset (in K)."
    bmetal_help = """Range of metallicities to include in the dataset. Half of the metallicity will
            C12, and the rest will be divided evenly among the remaining nuclei."""
    library_help = "Name of the library to use to load the network."
    targets_help = """Target nuclei to use for the algorithm. Will be protons and the endpoint by
            default."""
    tol_help = """Tolerance(s) to use for the algorithm. Will be [1e-3] + [1e-2]*(len(targets)-1)
            by default for 'drgep', and 1e-1 by default for 'pfa'."""
    no_mpi_help = "Disable MPI for this run."
    use_numpy_help = """If the algorithm has a 'use_numpy' option, turn it on. Some algorithms always
            run as if use_numpy=True."""
    first_pass_help = "If there is a first pass reduction option, turn it on."
    sa_help = "Error threshold to use when performing sensitivity analysis."
    
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-a', '--algorithm', default='drgep', help=algorithm_help)
    parser.add_argument('-e', '--endpoint', default=Nucleus('te108'), type=Nucleus,
            help=endpoint_help)
    parser.add_argument('-d', '--datadim', nargs=3, default=[4]*3, help=datadim_help)
    parser.add_argument('-dr', '--drho', type=int, help=drho_help)
    parser.add_argument('-dt', '--dtemp', help=dtemp_help)
    parser.add_argument('-dz', '--dmetal', type=int, help=dmetal_help)
    parser.add_argument('-br', '--brho', nargs=2, type=float, help=brho_help)
    parser.add_argument('-bt', '--btemp', nargs=2, type=float, help=btemp_help)
    parser.add_argument('-bz', '--bmetal', nargs=2, type=float, help=bmetal_help)
    parser.add_argument('-l', '--library', default="rp-process-lib", help=library_help)
    parser.add_argument('-t', '--targets', nargs='*', type=Nucleus, help=targets_help)
    parser.add_argument('--tol', nargs='*', type=float, help=tol_help)
    parser.add_argument('--no_mpi', action='store_true', help=no_mpi_help)
    parser.add_argument('--use_numpy', action='store_true', help=use_numpy_help)
    parser.add_argument('--first_pass', action='store_true', help=first_pass_help)
    parser.add_argument('-s', '--sens_analysis', default=0.05, type=float, help=sa_help)
    args = parser.parse_args(sys.argv[1:])
    
    args.algorithm = args.algorithm.lower()
    
    if args.algorithm == 'drgep':
        alg = drgep
        permute = (not args.use_numpy)
    elif args.algorithm == 'pfa':
        alg = pfa
        permute = False
    else:
        raise ValueError(f"Invalid algorithm selection: '{args.algorithm}'")
        
    if args.drho:
        args.datadim[0] = args.drho
    if args.dtemp:
        args.datadim[1] = args.dtemp
    if args.dmetal:
        args.datadim[2] = args.dmetal
        
    if not args.targets:
        args.targets = [Nucleus('p'), args.endpoint]
        
    if not args.tol:
        if args.algorithm == 'drgep':
            args.tol = [1e-3] + [1e-2] * (len(args.targets)-1)
        elif args.algorithm == 'pfa':
            args.tol = 0.9
    else:
        if args.algorithm == 'drgep':
            if len(args.tol) == 1:
                args.tol = [args.tol] * len(args.targets)
            elif len(args.tol) != len(args.targets):
                raise ValueError(f"For '{args.algorithm}', there should be one tolerance for each"
                        + " target nucleus.")
        elif args.algorithm == 'pfa':
            if len(args.tol) != 1:
                raise ValueError(f"For '{args.algorithm}', there should only be one tolerance.")
            else:
                args.tol = args.tol[0]
                
    #-----------------------------
    # Load the network and dataset
    #-----------------------------
    
    net = load_network(args.endpoint, library=args.library)
    data = list(dataset(net, args.datadim, permute, args.brho, args.btemp, args.bmetal))
    
    #-----------------------------
    # Prepare to run the algorithm
    #-----------------------------
    
    if args.algorithm == 'drgep':
        alg_args = \
        {
            'net': net,
            'conds': data,
            'targets': args.targets,
            'tols': args.tol,
            'returnobj': 'nuclei',
            'use_mpi': (not args.no_mpi),
            'use_numpy': args.use_numpy
        }
    elif args.algorithm == 'pfa':
        alg_args = \
        {
            'net': net,
            'conds': data,
            'targets': args.targets,
            'tol': args.tol,
            'returnobj': 'nuclei',
            'use_mpi': (not args.no_mpi),
            'first_pass': args.first_pass
        }
    
    if args.algorithm in {'pfa', 'drgep'}:
        args.algorithm = args.algorithm.upper()
    
    #------------------
    # Run the algorithm
    #------------------
    
    t0 = time.time()
    nuclei = alg(**alg_args)
    dt = time.time() - t0
    
    red_net = net.linking_nuclei(nuclei)
    second_data = list(dataset(net, args.datadim, True, args.brho, args.btemp, args.bmetal))
    errfunc = get_errfunc_enuc(net, second_data)
    
    if MPI.COMM_WORLD.Get_rank() == 0:
        
        print()
        print(f"{args.algorithm} reduction took {dt:.3f} s for {MPI.COMM_WORLD.Get_size()} processes.")
        print("Number of species in full network: ", len(net.unique_nuclei))
        print(f"Number of species in {args.algorithm} reduced network: ", len(nuclei))
        print("Reduced Network Error:", f"{errfunc(red_net)*100:.2f}%")
        print()
    
    # Perform sensitivity analysis
    MPI.COMM_WORLD.Barrier()
    t0 = time.time()
    # red_net = binary_search_trim(red_net, nuclei, errfunc)
    red_net = sens_analysis(red_net, errfunc, args.sens_analysis, True)
    dt = time.time() - t0
    
    if MPI.COMM_WORLD.Get_rank() == 0:
    
        print(f"SA reduction took {dt:.3f} s.")
        print(f"Number of species in {args.algorithm} + SA reduced network: ", len(red_net.unique_nuclei))
        print("Error: ", f"{errfunc(red_net)*100:.2f}%")
