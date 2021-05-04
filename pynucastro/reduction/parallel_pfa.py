import cProfile
import networkx as nx
import numpy as np
import pynucastro as pync
import sys
from collections import OrderedDict
from generate_data import dataset as get_conditions
from load_network import load_network
from mpi4py import MPI
from path_flux_analysis import calc_adj_matrix, graph_from_adj_matrix, get_remove_list, get_reduced_network
from pynucastro.networks import PythonNetwork
from pynucastro.rates import Library, RateFilter, Nucleus

def main(endpoint, targets =[Nucleus("p")], n=5, tol=0.4):
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    N_proc = comm.Get_size()

    net = load_network(endpoint)
    print("Network loaded on process %i" % rank)
    sys.stdout.flush()

    conds = list(get_conditions(net, n))
    if isinstance(n, int):
        n = np.ones(3, dtype=np.int32) * n
    else:
        n = np.array(list(n), dtype=np.int32)

    first_A = True
    n_conds = np.prod(n)
    for i in range(n_conds):
        if i % N_proc == rank:
            rho, T, comp = conds[i]
            rvals = net.evaluate_rates(rho=rho, T=T, composition=comp)
            if(not(i % (n_conds//10))):
                print("Proc %i on condition %i of %i" % (rank, i, n_conds))
                sys.stdout.flush()

            # grab adjacency matrix through PFA calculation on 2-neighbor paths
            A = calc_adj_matrix(net, rvals, tol)
            if first_A:
                A_red = np.copy(A)
                first_A = False
            else:
                A_red = np.maximum(A, A_red)

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
        G_pfa = graph_from_adj_matrix(net, A_final)
        r_species = get_remove_list(G_pfa, targets) # when working with many reaction conditions, intersection should be performed over all conditions

        # construct new network with only reactions involving reachable species
        reduced_net = get_reduced_network(net, r_species)

        print("Number of species in full network: ", len(net.unique_nuclei))
        print("Number of rates in full network: ", len(net.rates))
        print("Number of species in reduced network: ", len(reduced_net.unique_nuclei))
        print("Number of rates in reduced network: ", len(reduced_net.rates))

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
