import cProfile
import networkx as nx
import numpy as np
import pynucastro as pync
import sys
from collections import OrderedDict
from generate_data import dataset as get_conditions
from load_network import load_network
from mpi4py import MPI
import path_flux_analysis as pfa
from pynucastro.networks import PythonNetwork
from pynucastro.rates import Library, RateFilter, Nucleus

def main(endpoint, targets =[Nucleus("p")], n=5, tol=0.4):

    # Grab MPI settings and load data, conditions
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    N_proc = comm.Get_size()

    net = load_network(endpoint)
    print("Network loaded on process %i" % rank)
    sys.stdout.flush()

    conds = get_conditions(net, n)
    if isinstance(n, int):
        n = np.ones(3, dtype=np.int32) * n
    else:
        n = np.array(list(n), dtype=np.int32)

    rho_L = conds[0]
    T_L = conds[1]
    comp_L = conds[2]
    first_A = True
    n_conds = np.prod(n)

    # Precalculate data structures used in common over conditions
    r_map = pfa.get_r_map(net)
    r_set_indices = pfa.get_set_indices(net, r_map)
    stoich = pfa.get_stoich_matrix(net, r_map)

    # Iterate through conditions and reduce local matrix
    for i, rho in enumerate(rho_L):
        for j, T in enumerate(T_L):
            for k,comp in enumerate(comp_L):
                current = i*n[1]*n[2] + j*n[2] + k 
                if current % N_proc == rank:
                    rvals = net.evaluate_rates(rho=rho, T=T, composition=comp)
                    if(not(current % (n_conds//10))):
                        print("Proc %i on condition %i of %i" % (rank, current, n_conds))
                        sys.stdout.flush()

                    # grab adjacency matrix through PFA calculation on 2-neighbor paths
                    A = pfa.calc_adj_matrix(net, r_map, r_set_indices, stoich, rvals, tol)
                    if first_A:
                        A_red = np.copy(A)
                        first_A = False
                    else:
                        A_red = np.maximum(A, A_red)

    # Reduce adj matrix over processes
    A_final = np.zeros(A.shape)

    comm.Barrier()
    if(rank == 0):
        print("Reducing adjacency matrix across processes")
        sys.stdout.flush()
    comm.Reduce([A_red, MPI.DOUBLE], [A_final, MPI.DOUBLE], op=MPI.MAX, root=0)

    ## on rank 0 only
    # create new directed graph and perform DFS on targets to get nodes to remove
    if(rank==0):
        print("Reducing network.")
        sys.stdout.flush()
        G_pfa = pfa.graph_from_adj_matrix(net, A_final)
        r_species = pfa.get_remove_list(G_pfa, targets) # when working with many reaction conditions, intersection should be performed over all conditions

        # construct new network with only reactions involving reachable species
        reduced_net = pfa.get_reduced_network(net, r_species)

        print("Number of species in full network: ", len(net.unique_nuclei))
        print("Number of rates in full network: ", len(net.rates))
        print("Number of species in reduced network: ", len(reduced_net.unique_nuclei))
        print("Number of rates in reduced network: ", len(reduced_net.rates))
        print(reduced_net.unique_nuclei)

if __name__ == "__main__":
    if len(sys.argv) == 2:
        endpoint = Nucleus(sys.argv[1])
    elif len(sys.argv) == 1:
        endpoint = Nucleus('te108')
    else:
        print("Usage: ./load_network.py <endpoint>")

    pr = cProfile.Profile()
    pr.enable()
    main(endpoint)
    pr.disable()
    # Dump results:
    # - for binary dump
    comm = MPI.COMM_WORLD
    pr.dump_stats('cpu_%d.prof' % comm.rank)
    # - for text dump
    with open( 'cpu_%d.txt' % comm.rank, 'w') as output_file:
        sys.stdout = output_file
        pr.print_stats( sort='time' )
        sys.stdout = sys.__stdout__
