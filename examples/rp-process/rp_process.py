#!/usr/bin/env python3

import sys
import argparse
from collections import deque

from pynucastro.rates import Library, RateFilter
from pynucastro.nucdata import Nucleus, BindingTable
from pynucastro.networks import PythonNetwork, Composition

#################################################
#  Set up argument parser and process arguments #
#################################################
description = """Script for generating a reaction network modeling the rp-process with a given
        nucleus as the endpoint."""
        
endpoint_help = """The nucleus at which the network terminates. Should be provided as the string
        abbreviation of the nuclide (e.g. Ni56, case insensitive)."""
library_help = """The library file to draw the rates from. This is supplied directly to the Library
        constructor."""
write_network_help = """The name of the Python file to write the network to. No network will be
        written out if this is not supplied."""
write_lib_help = """The name of the file in the library directory to write the constructed library
        to. The library will not be written out if this is not supplied."""

parser = argparse.ArgumentParser(description=description)
parser.add_argument('endpoint', help=endpoint_help)
parser.add_argument('-l', '--library', default='reaclib_default2_20220329', help=library_help)
parser.add_argument('-w', '--write_network', default='', help=write_network_help)
parser.add_argument('--write_lib', default='', help=write_lib_help)
args = parser.parse_args(sys.argv[1:])

endpoint = Nucleus(args.endpoint)

#####################################
# Load library and generate network #
#####################################

print("Loading library...")

full_lib = Library(args.library)

print("Building network...")

core_nuclei = ["p", "d", "he3", "he4", "li7", "be7", "be8", "b8", "c12",
               "n13", "n14", "n15", "o14", "o15", "o16", "o17", "o18",
               "f17", "f18", "f19", "f20", "ne18", "ne19", "ne20", "ne21"]
core_nuclei = list(map(Nucleus, core_nuclei))
core_lib = full_lib.linking_nuclei(core_nuclei)

def is_beta_plus(rate):
    """ Filter for beta+ decays (and electron captures). """
    
    if len(rate.products) != len(rate.reactants):
        return False
    if len(rate.reactants) != 1:
        return False
    
    react, = rate.reactants
    prod, = rate.products
    return prod.Z < react.Z
    
# Restrict the library to these 7 reaction types

# Forward rates
p_gamma = RateFilter(reactants="p", max_products=1, exact=False)
alpha_gamma = RateFilter(reactants="he4", max_products=1, exact=False)
alpha_p = RateFilter(reactants="he4", products="p", exact=False)

# Reverse and weak rates
gamma_p = RateFilter(products="p", max_reactants=1, exact=False)
gamma_alpha = RateFilter(products="he4", max_reactants=1, exact=False)
p_alpha = RateFilter(reactants="p", products="he4", exact=False)
beta_plus = RateFilter(filter_function=is_beta_plus)

# Compute reduced library
red_lib = full_lib.filter((p_gamma, alpha_gamma, alpha_p,
        gamma_p, gamma_alpha, p_alpha, beta_plus))
        
bintable = BindingTable()

def flatten(iterable):
    """ Take iterable of iterables, and flatten it to one dimension. """
    
    for col in iterable:
        
        for item in col:
            
            yield item
            
def append_all(q, iterable):
    """ Append all items in the iterable to the queue. """
    
    for item in iterable: q.append(item)
        
def product_limiter():
    """
    This helps trim the library a bit by excluding rates with
    products with more protons than the endpoint, heavier than the
    endpoint, or with relatively high or low neutron percentages. 
    """
    
    # Proton number bounds
    Zlo, Zhi = 6, endpoint.Z
    # Nucleon number bounds
    Alo, Ahi = 12, endpoint.A
    # Bounds on A / Z ratio to drop peripheral nuclei
    Rlo, Rhi = 1.69, 2.2
    
    def limit_products(r):
        
        meet_conds = \
        (
            (Zlo <= p.Z <= Zhi and
            Alo <= p.A <= Ahi and
            Rlo <= p.A / p.Z <= Rhi and
            (p.N, p.Z) in bintable.energies) or
            (p.Z, p.A) == (1, 1) or
            (p.Z, p.A) == (2, 4)
            for p in r.products
        )
        return all(meet_conds)
    
    return limit_products
    
limiter = product_limiter()
final_lib = Library(rates=core_lib.get_rates())

seeds = [nuc for nuc in core_nuclei if nuc.A >= 12]
encountered = set(seeds) | {Nucleus("p"), Nucleus("he4")}
seeds = deque(seeds)

while seeds:
    
    # Get the new rates with seed as a reactant
    seed = seeds.popleft()
    filt = RateFilter(reactants=seed, filter_function=limiter, exact=False)
    new_lib = red_lib.filter(filt)
    if new_lib is None: continue
    final_lib += new_lib
    
    # Append all unseen nuclei to the queue
    prod = (r.products for r in new_lib.get_rates())
    prod = flatten(prod)
    prod = filter(lambda p: p not in encountered, prod)
    prod = sorted(set(prod))
    append_all(seeds, prod)
    encountered.update(prod)
    
encountered = sorted(encountered)
rp_net = PythonNetwork(libraries=[final_lib], precedence=["wc17", "ths8"])

print("Network constructed.")
print()
print(f"Species Encountered:")
print(encountered)
print()
print(f"Number of Species: {len(encountered)}")
print(f"Number of Rates: {len(rp_net.rates)}")
print()

if args.write_lib:
    final_lib.write_to_file(args.write_lib, True)

if args.write_network:
    print("Writing network...")
    rp_net.write_network(args.write_network)

print("Task completed.")
