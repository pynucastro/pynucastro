import numpy as np

from pynucastro import Composition, Nucleus

def calc_adj_matrix(net, rvals):
    
    N_species = len(net.unique_nuclei)
    
    # create index mapping
    j_map = dict()
    for idx, n in enumerate(net.unique_nuclei):
        j_map[n] = idx
    
    # Calculate coefficients for A
    p_A = np.zeros(N_species, dtype=np.float64)
    c_A = np.zeros(N_species, dtype=np.float64)
    
    # A along rows, B along columns
    p_AB = np.zeros((N_species, N_species))
    c_AB = np.zeros((N_species, N_species))
    
    for i in range(N_species):

        n = net.unique_nuclei[i]
        
        for r in net.nuclei_produced[n]:
            
            rval = r.products.count(n) * rvals[r]
            p_A[i] += rval
            
            bs = set(r.products) | set(r.reactants)
            for b in bs: p_AB[i,j_map[b]] += rval

        for r in net.nuclei_consumed[n]:
            
            rval = r.reactants.count(n) * rvals[r]
            c_A[i] += rval
            
            bs = set(r.products) | set(r.reactants)
            for b in bs: c_AB[i,j_map[b]] += rval

    denom = np.maximum(p_A, c_A)[:, np.newaxis]
    
    # Calculate direct interaction coefficients
    r_AB = np.abs(p_AB - c_AB) / denom

    return r_AB
    
def drgep_dijkstras(net, r_AB, target):
    
    # Number of species
    nspec = len(net.unique_nuclei)
    
    # Create data structures
    j_map = dict()
    
    inf = float('inf')
    
    dist = np.zeros(nspec, dtype=np.float64)
    R_TB = np.zeros(nspec, dtype=np.float64)
    
    for idx, n in enumerate(net.unique_nuclei):
        j_map[n] = idx
        dist[idx] = -inf
        R_TB[idx] = -inf
        
    dist[j_map[target]] = 1.0
    R_TB[j_map[target]] = inf
    
    imax = np.argmax(dist)
    
    # Main loop
    while dist[imax] > -inf:
        
        r_TA = dist[imax]
        R_TB[imax] = r_TA
        n = net.unique_nuclei[imax]
        bs = set()
        
        for r in net.nuclei_produced[n]:
            bs |= set(r.products) | set(r.reactants)
        for r in net.nuclei_consumed[n]:
            bs |= set(r.products) | set(r.reactants)
            
        for b in bs:
            jb = j_map[b]
            if R_TB[jb] > -inf:
                continue
            dist[jb] = max(dist[jb], r_TA * r_AB[imax, jb])
        
        dist[imax] = -inf
        imax = np.argmax(dist)
        
    return R_TB
    
def drgep(net, rvals, targets, tols):
    
    r_AB = calc_adj_matrix(net, rvals)
    R_TB = np.zeros_like(net.unique_nuclei, dtype=np.float64)
    
    try:
        targets = list(targets)
    except TypeError:
        targets = [targets]
        
    try:
        tols = list(tols)
    except TypeError:
        tols = [tols] * len(targets)
    
    for target, tol in zip(targets, tols):
        R_TB_i = drgep_dijkstras(net, r_AB, target)
        R_TB_i *= R_TB_i >= tol
        R_TB = np.maximum(R_TB, R_TB_i)

    nuclei = [net.unique_nuclei[i] for i in range(len(net.unique_nuclei))
              if R_TB[i] > 0.0]
    return net.linking_nuclei(nuclei)
    
if __name__ == "__main__":
    
    from pynucastro.reduction.load_network import load_network
    net = load_network(Nucleus('ni56'))
    
    T = 1.e9; rho = 1.e6
    comp = Composition(net.get_nuclei())
    comp.set_solar_like()
    rvals = net.evaluate_rates(rho=rho, T=T, composition=comp)

    targets = map(Nucleus, ['p', 'ni56'])
    reduced_net = drgep(net, rvals, targets, [1e-3, 1e-2])
    
    print("Number of species in full network: ", len(net.unique_nuclei))
    print("Number of rates in full network: ", len(net.rates))
    print("Number of species in reduced network: ", len(reduced_net.unique_nuclei))
    print("Number of rates in reduced network: ", len(reduced_net.rates))
    
    print(reduced_net.unique_nuclei)
