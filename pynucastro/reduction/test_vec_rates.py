import networkx as nx
import numpy as np
import pynucastro as pync
import sys
from load_network import load_network
from pynucastro.networks import PythonNetwork
from pynucastro.rates import Library, RateFilter, Nucleus
from collections import OrderedDict

def get_maps(net):
    n_map = dict()
    for i, n in enumerate(net.unique_nuclei):
        n_map[n] = i

    r_map = dict()
    for i, r in enumerate(net.rates):
        r_map[r] = i

    return n_map, r_map

def get_stoich_matrices(net, r_map):
    """
    Returns 3 matrices
    s_p, s_c, s_a
    s_p - coefs for rxns produced in, N_r x N_s
    s_c - coefs for rxns consumed in, N_r x N_s
    s_a - identities for rxns involved in, N-s x N_r
    """
    N_species = len(net.unique_nuclei)
    N_rates = len(net.rates)
    s_p = np.zeros((N_species, N_rates))
    s_c = np.zeros((N_species, N_rates))

    for i, n in enumerate(net.unique_nuclei):
        for r in net.nuclei_produced[n]:
            s_p[i, r_map[r]] = r.products.count(n)

        for r in net.nuclei_consumed[n]:
            s_c[i, r_map[r]] = r.reactants.count(n)

    return s_p, s_c, np.logical_or(s_p, s_c).astype(int).T


def main(endpoint, targets = [Nucleus("p")], tol=0.2):
    net = load_network(endpoint)
    print("Network loaded")

    T = 1e9
    rho = 1e4
    comp = pync.Composition(net.get_nuclei())
    comp.set_solar_like()
    pf_ref, yf_ref, rv_ref = net.evaluate_rates(rho=rho, T=T, composition=comp)
    pf_ref = np.array(pf_ref)
    yf_ref = np.array(yf_ref)
    rv_ref = np.array(rv_ref)

    n_map, r_map = get_maps(net)
    s_p, s_c, s_a = get_stoich_matrices(net, r_map)
    pf_test, yf_test, rv_test = net.evaluate_rates_arr(rho=rho, T=T, composition=comp, s_c=s_c)

    print(np.allclose(rv_ref, rv_test))
    print(np.allclose(pf_ref, pf_test))
    print(np.allclose(yf_ref, yf_test))


if __name__ == "__main__":
    if len(sys.argv) == 2:
        endpoint = Nucleus(sys.argv[1])
    elif len(sys.argv) == 1:
        endpoint = Nucleus('te108')
    else:
        print("Usage: ./load_network.py <endpoint>")
    main(endpoint)
