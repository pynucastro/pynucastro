#!/usr/bin/env python3

import sys

from pynucastro.networks import NumpyNetwork
from pynucastro.nucdata import Nucleus
from pynucastro.rates import Library, RateFilter


def load_network(endpoint=Nucleus('te108'), library_name='rp-process-lib'):
    """
    Load a network from a library, filtering the library so only nuclei with Z and A less than
    the endpoint's Z and A are included.
    """

    endpoint = Nucleus.cast(endpoint)

    def ff(rate):

        react_meet_conds = \
        (
            r.Z <= endpoint.Z and
            r.A <= endpoint.A
            for r in rate.reactants
        )

        prod_meet_conds = \
        (
            p.Z <= endpoint.Z and
            p.A <= endpoint.A
            for p in rate.products
        )

        return all(react_meet_conds) and all(prod_meet_conds)

    filt = RateFilter(filter_function=ff)
    lib = Library(library_name)
    lib = lib.filter(filt)
    return NumpyNetwork(libraries=[lib])


if __name__ == "__main__":

    args = {}
    if len(sys.argv) == 3:
        args = {"endpoint": Nucleus(sys.argv[1]), "library_name": sys.argv[2]}
    elif len(sys.argv) == 2:
        args = {"endpoint": Nucleus(sys.argv[1])}
    elif len(sys.argv) != 1:
        print("Usage: ./load_network.py [endpoint] [library_name]")
        sys.exit(0)

    net = load_network(**args)
    print("Network loaded!")
    print(f"Number of species: {len(net.unique_nuclei)}")
    print(f"Number of rates: {len(net.rates)}")
