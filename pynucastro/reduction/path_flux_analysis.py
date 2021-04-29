import networkx as nx
import numpy as np
import pynucastro as pync
import sys
from load_network import load_network
from pynucastro.networks import PythonNetwork
from pynucastro.rates import Library, RateFilter, Nucleus
from collections import OrderedDict

def first_pass_reduction(G, target_sources):

    # at this point, all the weights are 1
    paths = nx.multi_source_dijkstra_path(G, sources=target_sources, cutoff=2)
    u = []
    for p in paths:
        u = set(list(u) + list(paths[p]))

    # add only nodes accessible by distance of 2 from sources
    G_reduced = nx.DiGraph()
    G_reduced.add_nodes_from(u)

    edges = [e for e in G.edges if ((e[0] in u) and (e[1] in u))]
    G_reduced.add_edges_from(edges)

    return G_reduced

def calc_adj_matrix(net, rvals, tol):
    N_species = len(net.unique_nuclei)
    p_A = np.zeros(N_species)
    c_A = np.zeros(N_species)

    for i in range(N_species):
        #TODO: account for multiple nuclei being consumed/produced
        # example: he4 + he4 + he4 -> c12
        # example: he4 + he4 -> p + li7
        n = net.unique_nuclei[i]
        for r in net.nuclei_produced[n]:
            p_A[i] += rvals[r]

        for r in net.nuclei_consumed[n]:
            c_A[i] += rvals[r]

    denom = np.maximum(p_A, c_A)[:,np.newaxis]

    # create index mapping
    j_map = OrderedDict()
    nuclei_array = np.array(net.unique_nuclei)
    for n in net.unique_nuclei:
        j_map[n] = np.where(nuclei_array==n)[0][0]

    # A along rows, B along columns
    p_AB = np.zeros((N_species, N_species))
    c_AB = np.zeros((N_species, N_species))
    for i in range(N_species):
        n = net.unique_nuclei[i]
        for r in net.nuclei_produced[n]:
            for b in r.products:
                p_AB[i,j_map[b]] += rvals[r]

            for b in r.reactants:
                p_AB[i,j_map[b]] += rvals[r]

        for r in net.nuclei_consumed[n]:
            for b in r.products:
                c_AB[i,j_map[b]] += rvals[r]

            for b in r.reactants:
                c_AB[i,j_map[b]] += rvals[r]
    
    #by this point, should be in same form as pymars arrays
    r_pro_AB1 = p_AB/denom
    r_con_AB1 = c_AB/denom

    # direct copy of pymars PFA code with variable names changed
    r_pro_AB2 = np.zeros((N_species, N_species))
    r_con_AB2 = np.zeros((N_species, N_species))
    for i in range(N_species):
        pro1 = r_pro_AB1[:, i]
        pro2 = r_pro_AB1[i, :]
        con1 = r_con_AB1[:, i]
        con2 = r_con_AB1[i, :]
        pro1[i] = 0
        pro2[i] = 0
        con1[i] = 0
        con2[i] = 0
        r_pro_AB2 += np.outer(pro1, pro2)
        r_con_AB2 += np.outer(con1, con2)
    
    adjacency_matrix = r_pro_AB1 + r_con_AB1 + r_pro_AB2 + r_con_AB2
    np.fill_diagonal(adjacency_matrix, 0.0)

    # remove edges with weak dependence
    adjacency_matrix *= adjacency_matrix > tol
    return adjacency_matrix

def graph_from_adj_matrix(network, A):
    G = nx.from_numpy_matrix(A, create_using=nx.DiGraph())
    N_species = len(network.unique_nuclei)
    node_mapping = dict(zip([x for x in range(N_species)], network.unique_nuclei))

    return nx.relabel_nodes(G, node_mapping)

def get_remove_list(G, targets):
    lengths = nx.multi_source_dijkstra_path_length(G, targets)
    node_set = set(G.nodes)
    reachable = set(lengths.keys())
    return node_set.difference(reachable)

def get_reduced_network(net, r_species):
    r_rates = []
    for n in r_species:
        r_rates = set(list(r_rates) + list(net.nuclei_consumed[n]) + list(net.nuclei_produced[n]))
    
    return PythonNetwork(rates=list(set(net.rates).difference(r_rates)))

def main(endpoint, targets = [Nucleus("p")], tol=0.2):
    net = load_network(endpoint)
    print("Network loaded")
    # skipping the first pass reduction because it is mostly ineffectual
    # due to connectivity of p/he4
    T = 1e9
    rho = 1e4
    comp = pync.Composition(net.get_nuclei())
    comp.set_solar_like()
    rvals = net.evaluate_rates(rho=rho, T=T, composition=comp)

    # grab adjacency matrix through PFA calculation on 2-neighbor paths
    A = calc_adj_matrix(net, rvals, tol)

    # create new directed graph and perform DFS on targets to get nodes to remove
    G_pfa = graph_from_adj_matrix(net, A)
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
