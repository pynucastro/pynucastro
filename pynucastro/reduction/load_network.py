#!/usr/bin/env python3

import sys

from pynucastro.rates import Library, RateFilter
from pynucastro.nucdata import Nucleus
from pynucastro.networks import PythonNetwork

def load_network(endpoint=Nucleus('te108'), library='rp-process-lib'):
    
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
    lib = Library(library)
    lib = lib.filter(filt)
    return PythonNetwork(libraries=[lib])
    
if __name__ == "__main__":
    
    if len(sys.argv) == 2:
        endpoint = Nucleus(sys.argv[1])
    elif len(sys.argv) == 1:
        endpoint = Nucleus('te108')
    else:
        print("Usage: ./load_network.py <endpoint>")
        
    net = load_network(endpoint)
    print("Network loaded!")
    print(f"Number of species: {len(net.unique_nuclei)}")
    print(f"Number of rates: {len(net.rates)}")
