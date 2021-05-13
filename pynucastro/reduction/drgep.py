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
    
def drgep(net, conds, targets, tols, returnobj='net'):
    
    comm = MPI.COMM_WORLD
    MPI_N = comm.Get_size()
    MPI_rank = comm.Get_rank()
    
    # conds = _wrap_conds(conds)
    targets = _to_list(targets)
    tols = _to_list(tols, len(targets))
    
    R_TB_loc = np.zeros(len(net.unique_nuclei), dtype=np.float64)
    c_p, c_c, c_e = calc_count_matrices(net)

    n = np.array([len(conds[0]), len(conds[1]), len(conds[2])]).astype(int)
    rho_L = conds[0]
    T_L = conds[1]
    comp_L = conds[2]

    n_conds = np.prod(n)
    net.update_coef_arr()

    for k,comp in enumerate(comp_L):
        net.update_yfac_arr(composition=comp, s_c=c_c)

        for i, rho in enumerate(rho_L):
            net.update_prefac_arr(rho=rho, composition=comp)
            
            for j, T in enumerate(T_L):
                current = i*n[1]*n[2] + j*n[2] + k 
                if current % MPI_N == MPI_rank:
                    if(not(current % (n_conds//10))):
                        print("Proc %i on condition %i of %i" % (MPI_rank, current, n_conds))
                        sys.stdout.flush()

                    rvals_arr = net.evaluate_rates_arr(T=T)
                    # rvals = np.array(list(net.evaluate_rates(rho=rho, T=T, composition=comp).values()))
                    _drgep_kernel_numpy(net, R_TB_loc, c_p, c_c, c_e, rvals_arr, targets, tols)
        
    R_TB = np.zeros_like(R_TB_loc)
    comm.Allreduce([R_TB_loc, MPI.DOUBLE], [R_TB, MPI.DOUBLE], op=MPI.MAX)
              
    if returnobj == 'net':
        nuclei = [net.unique_nuclei[i] for i in range(len(net.unique_nuclei))
                  if R_TB[i] > 0.0]   
        return net.linking_nuclei(nuclei)
    elif returnobj == 'nuclei':
        idx = sorted(range(len(R_TB)), key=lambda i: R_TB[i], reverse=True)
        return [net.unique_nuclei[i] for i in idx if R_TB[i] > 0.0]
    else:
        raise ValueError("Invalid 'returnobj' argument: '{}'.".format(returnobj))
