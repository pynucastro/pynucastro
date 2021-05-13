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
import time


def main(endpoint, targets =[Nucleus("p")], n=16, tol=0.4):

    # Grab MPI settings and load data, conditions
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    N_proc = comm.Get_size()

    net = load_network(endpoint)
    # print("Network loaded on process %i" % rank)
    sys.stdout.flush()

    conds = get_conditions(net, n)
    if isinstance(n, int):
        n = np.ones(3, dtype=np.int32) * n
    else:
        n = np.array(list(n), dtype=np.int32)

    #only designing this to work for 2^n processes, 2^m conditions for m >= 4
    if(n[2] >= N_proc):
        comp_i = (n[2]//N_proc)*rank
        comp_f = (n[2]//N_proc)*(rank+1)
        rho_i = 0
        rho_f = -1
    else: 
        comp_i = rank % n[2] 
        comp_f = comp_i + 1
        comp_split = N_proc / n[2]
        rho_i = (n[0]//(N_proc / n[2])) * (rank// n[2])
        rho_f = (n[0]//(N_proc / n[2])) * (rank// n[2] + 1)


    rho_L = conds[0][rho_i:rho_f]
    T_L = conds[1]
    comp_L = conds[2][comp_i:comp_f]

    first_A = True
    n_conds = np.prod(n)

    comm.Barrier() #synchronize timing
    t_0 = MPI.Wtime()
    # Precalculate data structures used in common over conditions
    n_map, r_map = pfa.get_maps(net)
    r_set_indices = pfa.get_set_indices(net, n_map)
    s_p, s_c, s_a = pfa.get_stoich_matrices(net, r_map)

    net.update_coef_arr()

    t_1 = MPI.Wtime()
    # Iterate through conditions and reduce local matrix
    count = 0
    for k,comp in enumerate(comp_L):
        net.update_yfac_arr(composition=comp, s_c=s_c)
        for i, rho in enumerate(rho_L):
            net.update_prefac_arr(rho=rho, composition=comp)
            for j, T in enumerate(T_L):
                rvals_arr = net.evaluate_rates_arr(T=T)
                # if(not(current % (n_conds//10))):
                #     print("Proc %i on condition %i of %i" % (rank, current, n_conds))
                #     sys.stdout.flush()

                # grab adjacency matrix through PFA calculation on 2-neighbor paths
                A = pfa.calc_adj_matrix(net, s_p, s_c, s_a, rvals_arr)
                if first_A:
                    A_red = np.copy(A)
                    first_A = False
                else:
                    A_red = np.maximum(A, A_red)

    np.fill_diagonal(A_red, 0.0)

    # remove edges with weak dependence
    A_red *= A_red > tol

    # Reduce adj matrix over processes
    A_final = np.zeros(A.shape)


    t_2 = MPI.Wtime()
    comm.Barrier()
    t_B = MPI.Wtime()
    # if(rank == 0):
    #     print("Reducing adjacency matrix across processes")
    #     sys.stdout.flush()
    comm.Reduce([A_red, MPI.DOUBLE], [A_final, MPI.DOUBLE], op=MPI.MAX, root=0)
    t_3 = MPI.Wtime()

    ## on rank 0 only
    # create new directed graph and perform DFS on targets to get nodes to remove
    print("Process %i calculated %i conditions" % (rank, count))
    sys.stdout.flush()
    comm.Barrier()
    t_0_arr = np.array(comm.gather(t_0, root=0))
    t_1_arr = np.array(comm.gather(t_1, root=0))
    t_2_arr = np.array(comm.gather(t_2, root=0))
    t_B_arr = np.array(comm.gather(t_B, root=0))
    t_3_arr = np.array(comm.gather(t_3, root=0))
    if(rank==0):
        G_pfa = pfa.graph_from_adj_matrix(net, A_final)
        r_species = pfa.get_remove_list(G_pfa, targets) # when working with many reaction conditions, intersection should be performed over all conditions

        t_4 = MPI.Wtime()
        # construct new network with only reactions involving reachable species
        reduced_net = pfa.get_reduced_network(net, r_species)
        t_5 = MPI.Wtime()
        print("Number of species in full network: ", len(net.unique_nuclei))
        print("Number of rates in full network: ", len(net.rates))
        print("Number of species in reduced network: ", len(reduced_net.unique_nuclei))
        print("Number of rates in reduced network: ", len(reduced_net.rates))
        print(f"PFA took {t_5-t_0:.3f} s overall.")

        dt1 = t_1_arr-t_0_arr
        dt2 = t_2_arr-t_1_arr
        dt3 = t_B_arr-t_2_arr
        dt4 = t_3_arr-t_B_arr
        print(f"Precalculation mean: {np.mean(dt1):.3f} s std. dev:  {np.std(dt1):.3f}.")
        print(f"Local adj matrix mean: {np.mean(dt2):.3f} s std. dev:  {np.std(dt2):.3f}.")
        print(f"Barrier waiting mean: {np.mean(dt3):.3f} s std. dev:  {np.std(dt3):.3f}.")
        print(f"Parallel reduction mean: {np.mean(dt4):.3f} s std. dev:  {np.std(dt4):.3f}.")
        print(f"{t_4-t_3:.3f} s in graph traversal.")
        print(f"{t_5-t_4:.3f} s in forming final network.")
        # print(reduced_net.unique_nuclei)

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
    # pr.dump_stats('cpu_%d.prof' % comm.rank)
    # - for text dump
    with open( 'cpu_%d.txt' % comm.rank, 'w') as output_file:
        sys.stdout = output_file
        pr.print_stats( sort='time' )
        sys.stdout = sys.__stdout__
