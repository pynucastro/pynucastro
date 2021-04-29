import networkx as nx
import numpy as np
import pynucastro as pync
import sys
from load_network import load_network
from pynucastro.networks import PythonNetwork
from pynucastro.rates import Library, RateFilter, Nucleus
from collections import OrderedDict
from mpi4py import MPI
from path_flux_analysis import calc_adj_matrix, graph_from_adj_matrix, get_remove_list, get_reduced_network
from generate_data import dataset as get_conditions

def main(endpoint, targets =[Nucleus("p")], n=10, tol=0.2):
    net = load_network(endpoint)
    print("Network loaded")

    conds = get_conditions(net, n)
    if isinstance(n, int):
        n = np.ones(3, dtype=np.int32) * n
    else:
        n = np.array(list(n), dtype=np.int32)

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    N_proc = comm.Get_size()


    rho_L = conds[0]
    T_L = conds[1]
    comp_L = conds[2]
    first_A = True
    for i, rho in enumerate(rho_L):
        for j, T in enumerate(T_L):
            for k,comp in enumerate(comp_L):
                iter = i*n[1]*n[2] + j*n[2] + k 
                if iter % N_proc == rank:
                    rvals = net.evaluate_rates(rho=rho, T=T, composition=comp)

                    # grab adjacency matrix through PFA calculation on 2-neighbor paths
                    A = calc_adj_matrix(net, rvals, tol)
                    if first_A:
                        A_red = np.copy(A)
                        first_A = False
                    else:
                        A_red = np.maximum(A, A_red)

    if(rank ==0):
        A_final = np.zeros(A.shape)

    comm.Barrier()
    comm.Reduce([A_red, MPI.DOUBLE], [A_final, MPI.DOUBLE], op=MPI.MAX, root=0)

    ## on rank 0 only
    # create new directed graph and perform DFS on targets to get nodes to remove
    if(rank==0):
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
    main(endpoint)
