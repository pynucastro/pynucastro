import numpy as np

from collections import namedtuple
from pynucastro import Composition, Nucleus
from pynucastro.nucdata import BindingTable
from mpi4py import MPI
import sys

def calc_count_matrices(net):
    
    # Rate -> index mapping
    r_map = dict()
    for i, r in enumerate(net.rates):
        r_map[r] = i
    
    N_species = len(net.unique_nuclei)
    N_rates = len(net.rates)
    
    # Counts for reactions producing nucleus
    c_p = np.zeros((N_species, N_rates), dtype=np.int32)
    # Counts for reactions consuming nucleus
    c_c = np.zeros((N_species, N_rates), dtype=np.int32)

    for i, n in enumerate(net.unique_nuclei):
        
        for r in net.nuclei_produced[n]:
            c_p[i, r_map[r]] = r.products.count(n)

        for r in net.nuclei_consumed[n]:
            c_c[i, r_map[r]] = r.reactants.count(n)
    
    # Whether the nucleus is involved in the reaction or not
    c_e = np.logical_or(c_p, c_c).astype(np.int32).T
            
    return c_p, c_c, c_e

def calc_adj_matrix_numpy(net, c_p, c_c, c_e, rvals_arr):
    
    # Evaluate terms on RHS of ODE system
    prod_terms = c_p * rvals_arr
    cons_terms = c_c * rvals_arr
    
    # Calculate total production and consumption of each nucleus A
    p_A = prod_terms.sum(axis=1)
    c_A = cons_terms.sum(axis=1)
    
    # Calculate production / consumption of A in reactions involving B
    p_AB = prod_terms @ c_e
    c_AB = cons_terms @ c_e
    
    # We will normalize by maximum of production and consumption fluxes
    denom = np.maximum(p_A, c_A)[:, np.newaxis]
    
    # Calculate direct interaction coefficients
    r_AB = np.abs(p_AB - c_AB) / denom
    
    return r_AB
    
def calc_adj_matrix(net, rvals):
    
    N_species = len(net.unique_nuclei)
    
    # create index mapping
    j_map = dict()
    for idx, n in enumerate(net.unique_nuclei):
        j_map[n] = idx
    
    # Calculate coefficients for A
    p_A = np.zeros(N_species, dtype=np.float64)
    c_A = np.zeros(N_species, dtype=np.float64)
    
    # A along rows, B along columns
    p_AB = np.zeros((N_species, N_species))
    c_AB = np.zeros((N_species, N_species))
    
    for i in range(N_species):

        n = net.unique_nuclei[i]
        
        for r in net.nuclei_produced[n]:
            
            rval = r.products.count(n) * rvals[r]
            p_A[i] += rval
            
            bs = set(r.products) | set(r.reactants)
            for b in bs: p_AB[i,j_map[b]] += rval

        for r in net.nuclei_consumed[n]:
            
            rval = r.reactants.count(n) * rvals[r]
            c_A[i] += rval
            
            bs = set(r.products) | set(r.reactants)
            for b in bs: c_AB[i,j_map[b]] += rval

    denom = np.maximum(p_A, c_A)[:, np.newaxis]
    
    # Calculate direct interaction coefficients
    r_AB = np.abs(p_AB - c_AB) / denom

    return r_AB
    
def drgep_dijkstras(net, r_AB, target):
    
    # Number of species
    nspec = len(net.unique_nuclei)
    
    # Create data structures
    j_map = dict()
    
    inf = float('inf')
    
    dist = np.zeros(nspec, dtype=np.float64)
    R_TB = np.zeros(nspec, dtype=np.float64)
    
    for idx, n in enumerate(net.unique_nuclei):
        j_map[n] = idx
        dist[idx] = -inf
        R_TB[idx] = -inf
        
    dist[j_map[target]] = 1.0
    R_TB[j_map[target]] = inf
    
    imax = np.argmax(dist)
    
    # Main loop
    while dist[imax] > -inf:
        
        r_TA = dist[imax]
        R_TB[imax] = r_TA
        n = net.unique_nuclei[imax]
        bs = set()
        
        for r in net.nuclei_produced[n]:
            bs |= set(r.products) | set(r.reactants)
        for r in net.nuclei_consumed[n]:
            bs |= set(r.products) | set(r.reactants)
            
        for b in bs:
            jb = j_map[b]
            if R_TB[jb] > -inf:
                continue
            dist[jb] = max(dist[jb], r_TA * r_AB[imax, jb])
        
        dist[imax] = -inf
        imax = np.argmax(dist)
        
    return R_TB

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
    
def _to_list(x, n=1):
    
    try:
        return list(x)
    except TypeError:
        return [x] * n
    
def _drgep_kernel(net, R_TB, rvals, targets, tols):
    
    r_AB = calc_adj_matrix(net, rvals)
    
    for target, tol in zip(targets, tols):
        
        R_TB_i = drgep_dijkstras(net, r_AB, target)
        np.maximum(R_TB, R_TB_i, out=R_TB, where=(R_TB_i >= tol))

def _drgep_kernel_numpy(net, R_TB, c_p, c_c, c_e, rvals, targets, tols):
    
    r_AB = calc_adj_matrix_numpy(net, c_p, c_c, c_e, rvals)
    
    for target, tol in zip(targets, tols):
        
        R_TB_i = drgep_dijkstras(net, r_AB, target)
        np.maximum(R_TB, R_TB_i, out=R_TB, where=(R_TB_i >= tol))
    
def drgep(net, conds, targets, tols):
    
    comm = MPI.COMM_WORLD
    MPI_N = comm.Get_size()
    MPI_rank = comm.Get_rank()
    
    conds = _wrap_conds(conds)
    targets = _to_list(targets)
    tols = _to_list(tols, len(targets))
    
    R_TB_loc = np.zeros(len(net.unique_nuclei), dtype=np.float64)
    c_p, c_c, c_e = calc_count_matrices(net)
    
    for i in range(MPI_rank, len(conds), MPI_N):
        if(not(i % (len(conds)//10))):
            print("Proc %i on condition %i of %i" % (MPI_rank, i, len(conds)))
            sys.stdout.flush()

        rho, T, comp = conds[i]
        rvals = np.array(list(net.evaluate_rates(rho=rho, T=T, composition=comp).values()))
        _drgep_kernel_numpy(net, R_TB_loc, c_p, c_c, c_e, rvals, targets, tols)
        
    R_TB = np.zeros_like(R_TB_loc)
    comm.Allreduce([R_TB_loc, MPI.DOUBLE], [R_TB, MPI.DOUBLE], op=MPI.MAX)
        
    nuclei = [net.unique_nuclei[i] for i in range(len(net.unique_nuclei))
              if R_TB[i] > 0.0]
    return net.linking_nuclei(nuclei)

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
        ebind[i] = bintable[n.N, n.Z].nucbind
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

def rel_err(x, x0):

    return np.abs((x - x0) / x0)

def error_function(net_new, net_old, conds):

    conds = _wrap_conds(conds)
    err = np.zeros(3, dtype=np.float64)

    for rho, T, comp in conds:

        # comp_red = map_comp2(comp, net_new)

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
    
if __name__ == "__main__":
    
    from pynucastro.reduction.load_network import load_network
    from pynucastro.reduction.generate_data import dataset
    import time
    
    net = load_network(Nucleus('ni56'))
    data = list(dataset(net, n=4))
    
    targets = map(Nucleus, ['p', 'ni56'])
    t0 = time.time()
    reduced_net = drgep(net, data, targets, [1e-3, 1e-2])
    dt = time.time() - t0
    
    if MPI.COMM_WORLD.Get_rank() == 0:
        print(f"DRGEP reduction took {dt:.3f} s.")
        print("Number of species in full network: ", len(net.unique_nuclei))
        print("Number of rates in full network: ", len(net.rates))
        print("Number of species in reduced network: ", len(reduced_net.unique_nuclei))
        print("Number of rates in reduced network: ", len(reduced_net.rates))
        
        # print()
        # print("Evaluating error (max across all sets of conditions)...")
        # err = error_function(reduced_net, net, data)
        # print("enuc_dot: {:.1f}%; ye_dot: {:.1f}%; abar_dot: {:.1f}%".format(*(100*err)))
