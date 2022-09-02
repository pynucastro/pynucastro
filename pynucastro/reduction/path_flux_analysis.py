import numpy as np
import networkx as nx
from .reduction_utils import FailedMPIImport, mpi_numpy_decomp
try:
    from mpi4py import MPI
except ImportError:
    MPI = FailedMPIImport()

def pfa_first_pass_reduction(net, target_sources):
    """
    Takes a network and removes all nodes that are more than two edges away from the
    nuclei in *target_sources*.
    """

    # at this point, all the weights are 1
    G = net.get_networkx_graph()
    paths = nx.multi_source_dijkstra_path(G, sources=target_sources, cutoff=2)
    
    # add only nodes accessible by distance of 2 from sources
    u = set()
    for p in paths:
        u = set(list(u) + list(paths[p]))
    
    return net.linking_nuclei(u)

def calc_adj_matrix(net, rvals_arr):
    """Calculate adjacency graph where edge weights are PFA interaction coefficients."""
    
    s_p_scaled = net.nuc_prod_count*rvals_arr
    s_c_scaled = net.nuc_cons_count*rvals_arr

    p_A = np.sum(s_p_scaled, axis=1)
    c_A = np.sum(s_c_scaled, axis=1)
    denom = np.maximum(p_A, c_A)[:,np.newaxis]

    p_AB = s_p_scaled @ net.nuc_used
    c_AB = s_c_scaled @ net.nuc_used
    r_pro_AB1 = p_AB/denom
    r_con_AB1 = c_AB/denom

    np.fill_diagonal(r_pro_AB1, 0.0)
    np.fill_diagonal(r_con_AB1, 0.0)

    r_pro_AB2 = r_pro_AB1 @ r_pro_AB1
    r_con_AB2 = r_con_AB1 @ r_con_AB1
    
    adjacency_matrix = r_pro_AB1 + r_con_AB1 + r_pro_AB2 + r_con_AB2
    np.fill_diagonal(adjacency_matrix, 0.0)

    return adjacency_matrix

def graph_from_adj_matrix(network, A):
    """Convert PFA adjacency matrix to networkx graph."""
    
    G = nx.from_numpy_matrix(A, create_using=nx.DiGraph())
    N_species = len(network.unique_nuclei)
    node_mapping = dict(zip([x for x in range(N_species)], network.unique_nuclei))

    return nx.relabel_nodes(G, node_mapping)

def get_reachable_nodes(G, targets):
    """Get list of nodes (nuclei) reachable from a set of targets."""
    
    lengths = nx.multi_source_dijkstra_path_length(G, targets)
    reachable = list(lengths.keys())
    return reachable
    
def reduced_net_from_adj_matrix(basenet, adj, targets):
    """
    Calculate and return reduced network from original net, PFA adjacency matrix, and target nuclei.
    """
    
    G_pfa = graph_from_adj_matrix(basenet, adj)
    reachspec = get_reachable_nodes(G_pfa, targets)
    return basenet.linking_nuclei(reachspec)
    
def _pfa(net, conds):

    #----------------------------------------
    # Unpack conditions; make precalculations
    #----------------------------------------
    
    comp_L, rho_L, T_L = conds
    
    net.calc_count_matrices()
    net.update_rate_coef_arr()
    
    #------------------------------------------------
    # Calculate interaction coefficients (vectorized)
    #------------------------------------------------
    
    A = np.zeros((len(net.unique_nuclei), len(net.unique_nuclei)) dtype=np.float64)

    for comp in comp_L:
        net.update_yfac_arr(comp)
        for rho in rho_L:
            net.update_prefac_arr(rho, comp)
            for T in T_L:
                rvals_arr = net.evaluate_rates_arr(T)
                A = np.maximum(A, calc_adj_matrix(net, rvals_arr))
    
    net.clear_arrays()
    return A
    
def _pfa_mpi(net, conds):

    #----------------------------------------
    # Unpack conditions; make precalculations
    #----------------------------------------
    
    n = tuple(map(len, conds))
    comp_L, rho_L, T_L = conds
    
    net.calc_count_matrices()
    net.update_rate_coef_arr()
    
    #-----------------------------------
    # Init. MPI and divide up conditions
    #-----------------------------------
    
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    N_proc = comm.Get_size()
    
    A = np.zeros((len(net.unique_nuclei), len(net.unique_nuclei)) dtype=np.float64)
    comp_idx, comp_step, rho_idx, rho_step, T_idx, T_step = mpi_numpy_decomp(N_proc, rank, n)

    # Iterate through conditions and reduce local matrix
    for i in range(comp_idx, n[0], comp_step):
        comp = comp_L[i]
        net.update_yfac_arr(comp)
        for j in range(rho_idx, n[1], rho_step):
            rho = rho_L[j]
            net.update_prefac_arr(rho, comp)
            for k in range(T_idx, n[2], T_step):
                T = T_L[k]
                rvals_arr = net.evaluate_rates_arr(T)
                # grab adjacency matrix through PFA calculation on 2-neighbor paths
                A = np.maximum(A, calc_adj_matrix(net, rvals_arr))
                
    # Reduce adj matrix over processes
    A_final = np.zeros(A.shape)

    # if SA:
    #     comm.Allreduce([A_red, MPI.DOUBLE], [A_final, MPI.DOUBLE], op=MPI.MAX, root=0)
    # else:
    #     comm.Reduce([A_red, MPI.DOUBLE], [A_final, MPI.DOUBLE], op=MPI.MAX, root=0)
    comm.Allreduce([A, MPI.DOUBLE], [A_final, MPI.DOUBLE], op=MPI.MAX)
    
    net.clear_arrays()
    return A_final
    
def pfa(net, conds, targets=None, tol=None, returnobj='net', use_mpi=False, first_pass=False):
    """
    Implementation of Path Flux Analysis (PFA) reduction method described in Sun et al. 2010
    (doi:10.1016/j.combustflame.2010.03.006).
    
    :param net: The network (RateCollection) to reduce.
    :param conds: A set of conditions to reduce over. Should either be a sequence of (composition,
        density, temperature) sequences/tuples if running in standard mode, or a sequence of 3
        sequences ((composition, density, temperature) ordering) if running in NumPy mode. In the
        latter case, the sequences will be permuted to create the dataset. The compositions should
        be pynucastro Composition objects.
    :param targets: A collection of target nuclei (or a single target nucleus) to run the
        graph search algorithm from. Should be supplied as pynucastro Nucleus objects. Note that
        these are needed to get the reduced net, but do not impact any of the interaction
        coefficients.
    :param tol: Interaction coefficient tolerance to use when removing edges from the network. Only
        required if requesting 'net' or 'nuclei' as the *returnobj*. If coefficients are requested
        and this parameter is supplied, edge coefficents that are below the tolerance will be set
        to 0 before returning the edge coefficients or calculating the nucleus coefficients.
    :param returnobj: The type of object to return. Options are 'net' (a reduced network, the default
        setting), 'nuclei' (unique nuclei in the reduced network, ordered so the interaction
        coefficients are descending), 'edge_coeff' (the interaction coefficients for each edge
        connecting a pair of nuclei in the original net), and 'nuc_coeff' (the sum of edge
        interaction coefficients for each nucleus in the original net).
    :param use_mpi: Whether to divide up the set of conditions across MPI processes or not. Default
        setting is *False*.
    :param use_numpy: Whether to use NumPy to vectorize the interaction coefficient calculations or
        not. This is more memory intensive and may actually hinder performance for some setups.
        Conditions should be supplied as 3 lists that will be permuted to form the dataset (see
        *conds* parameter). Default setting is *False*.
    """
    
    #------------------
    # Process arguments
    #------------------
    
    if returnobj not in {'net', 'nuclei', 'edge_coeff', 'nuc_coeff'}:
        raise ValueError(f"Invalid 'returnobj' argument: '{returnobj}'.")
        
    can_reduce = (targets is not None) and (tol is not None)
    if returnobj == 'net' or returnobj == 'nuclei':
        assert can_reduce, "Must supply both 'targets' and 'tol' as arguments if requesting" +
                f" returnobj '{returnobj}'."
        
    targets = _to_list(targets)
    
    #-------------------------------------------------------
    # Determine operation mode and launch appropriate helper
    #-------------------------------------------------------
    
    if use_mpi:
        adj = _pfa_mpi(net, conds)
    else:
        adj = _pfa(net, conds)
        
    #------------------------
    # Return requested object
    #------------------------
    
    if tol is not None:
        adj *= adj > tol
    
    if returnobj == 'net':
        return reduced_net_from_adj_matrix(net, adj, targets)
    elif returnobj == 'nuclei':
        R_TB = np.sum(adj, axis=0)
        nuc = set(reduced_net_from_adj_matrix(net, adj, targets).get_nuclei())
        idx = sorted(range(len(R_TB)), key=lambda i: R_TB[i], reverse=True)
        return [net.unique_nuclei[i] for i in idx if net.unique_nuclei[i] in nuc]
    elif returnobj == 'edge_coeff':
        return adj
    elif returnobj == 'nuc_coeff':
        return np.sum(adj, axis=0)
