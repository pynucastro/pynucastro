import networkx as nx
from load_network import load_network
from pynucastro.rates import Library, RateFilter, Nucleus
from pynucastro.networks import PythonNetwork

def first_pass_reduction(G, target_sources):

    # at this point, all the weights are 1
    paths = nx.multi_source_dijkstra_path_length(G, sources=target_sources, cutoff=2)
    u = []
    for p in paths:
        u = set(list(u) + list(paths[p]))

    # add only nodes accessible by distance of 2 from sources
    G_reduced = nx.DiGraph()
    G_reduced.add_nodes_from(u)

    edges = [e for e in G.edges if ((e[0] in u) and (e[1] in u))]
    G_reduced.add_edges_from(edges)

    return G_reduced


def main(endpoint):
    net = load_network(endpoint)
    G = net.get_reaction_network_graph()
    print("Original order %i " % G.order())

    target_sources = ["p", "he4", "li7"]
    G_reduced = first_pass_reduction(G)
    print("Reduced order %i " % G_reduced.order())


if __name__ == "__main__":
    if len(sys.argv) == 2:
        endpoint = Nucleus(sys.argv[1])
    elif len(sys.argv) == 1:
        endpoint = Nucleus('te108')
    else:
        print("Usage: ./load_network.py <endpoint>")
    main(endpoint)