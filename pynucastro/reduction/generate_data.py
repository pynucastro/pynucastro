#!/usr/bin/env python3

from pynucastro import Composition, Nucleus
import numpy as np

def dataset(network, n=10):
    
    if isinstance(n, int):
        n = np.ones(3, dtype=np.int32) * n
    else:
        n = np.array(list(n), dtype=np.int32)
    
    # Bounds on each variable
    b_rho = (1e2, 1e6) # density (g/cm^3)
    b_T = (8e6, 1.5e9) # temperature (K)
    b_Z = (0.02, 0.2) # metallicity
    
    rho = np.logspace(*map(np.log10, b_rho), num=n[0])
    T = np.logspace(*map(np.log10, b_T), num=n[1])
    comp_list = []
    
    for Z in np.linspace(*b_Z, num=n[2]):
        
        comp = Composition(network.get_nuclei())
        # 75/25 H1/He4 ratio
        comp.X[Nucleus("p")] = (1. - Z) * 0.75
        comp.X[Nucleus("he4")] = (1. - Z) * 0.25
        # 50% of remaining stuff in carbon-12
        comp.X[Nucleus("c12")] = Z * 0.5
        # Split the rest across all other nuclei
        omit = {Nucleus("p"), Nucleus("he4"), Nucleus("c12")}
        rem = 0.5 * Z / (len(comp.X) - len(omit))
        
        for nuc in comp.X:
            if nuc not in omit:
                comp.X[nuc] = rem
                
        comp.normalize()
        comp_list.append(comp)
                
    for rho_i in rho:
        for T_i in T:
            for comp in comp_list:
                yield (rho_i, T_i, comp)

if __name__ == "__main__":
    
    from pynucastro.reduction.load_network import load_network
    network = load_network(Nucleus("ni56"))
    conds_list = list(dataset(network))
    rho = sorted({conds[0] for conds in conds_list})
    T = sorted({conds[1] for conds in conds_list})
    comp = {conds[2] for conds in conds_list}
    
    print("ρ")
    print(rho)
    print()
    
    print("T")
    print(T)
    print()
    
    def comp_converter(comp):
        
        X = comp.X[Nucleus("p")]
        Y = comp.X[Nucleus("he4")]
        ZC12 = comp.X[Nucleus("c12")]
        Zmin = min(comp.X.values())
        return (X, Y, ZC12, Zmin)
        
    print("X")
    print(np.array(list(map(comp_converter, comp))))
