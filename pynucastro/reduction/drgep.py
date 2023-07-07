import numpy as np

from pynucastro.reduction.reduction_utils import (mpi_importer,
                                                  mpi_numpy_decomp, to_list)

MPI = mpi_importer()


def calc_interaction_matrix(net, rvals):
    """Calculate direct interaction coefficients."""

    N_species = len(net.unique_nuclei)

    # create index mapping
    j_map = {}
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
            for b in bs:
                p_AB[i, j_map[b]] += rval

        for r in net.nuclei_consumed[n]:

            rval = r.reactants.count(n) * rvals[r]
            c_A[i] += rval

            bs = set(r.products) | set(r.reactants)
            for b in bs:
                c_AB[i, j_map[b]] += rval

    denom = np.maximum(p_A, c_A)[:, np.newaxis]

    # Calculate direct interaction coefficients
    r_AB = np.abs(p_AB - c_AB) / denom

    return r_AB


def calc_interaction_matrix_numpy(net, rvals_arr):
    """Calculate direct interaction coefficients using NumPy."""

    # Evaluate terms on RHS of ODE system
    prod_terms = net.nuc_prod_count * rvals_arr
    cons_terms = net.nuc_cons_count * rvals_arr

    # Calculate total production and consumption of each nucleus A
    p_A = prod_terms.sum(axis=1)
    c_A = cons_terms.sum(axis=1)

    # Calculate production / consumption of A in reactions involving B
    p_AB = prod_terms @ net.nuc_used
    c_AB = cons_terms @ net.nuc_used

    # We will normalize by maximum of production and consumption fluxes
    denom = np.maximum(p_A, c_A)[:, np.newaxis]

    # Calculate direct interaction coefficients
    r_AB = np.abs(p_AB - c_AB) / denom

    return r_AB


def get_adj_nuc(net):
    """
    Get set of adjacent nuclei for each nucleus in the net. Returns dictionary keyed by
    nuclei.
    """

    adj_nuc = {}

    for n in net.unique_nuclei:

        bs = set()

        for r in net.nuclei_produced[n]:
            bs |= set(r.products) | set(r.reactants)
        for r in net.nuclei_consumed[n]:
            bs |= set(r.products) | set(r.reactants)

        adj_nuc[n] = bs

    return adj_nuc


def drgep_dijkstras(net, r_AB, target, adj_nuc):
    """
    Modified Dijkstra's algorithm to find paths that maximize the overall interaction coefficient
    between the target and each other nucleus in the net.
    """

    # Number of species
    nspec = len(net.unique_nuclei)

    # Create data structures
    j_map = {}

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
        bs = adj_nuc[n]

        for b in bs:
            jb = j_map[b]
            if R_TB[jb] > -inf:
                continue
            dist[jb] = max(dist[jb], r_TA * r_AB[imax, jb])

        dist[imax] = -inf
        imax = np.argmax(dist)

    return R_TB


def _drgep_kernel(net, R_TB, rvals, targets, tols, adj_nuc):

    r_AB = calc_interaction_matrix(net, rvals)

    for target, tol in zip(targets, tols):

        R_TB_i = drgep_dijkstras(net, r_AB, target, adj_nuc)
        np.maximum(R_TB, R_TB_i, out=R_TB, where=R_TB_i >= tol)


def _drgep_kernel_numpy(net, R_TB, rvals, targets, tols, adj_nuc):

    r_AB = calc_interaction_matrix_numpy(net, rvals)

    for target, tol in zip(targets, tols):

        R_TB_i = drgep_dijkstras(net, r_AB, target, adj_nuc)
        np.maximum(R_TB, R_TB_i, out=R_TB, where=R_TB_i >= tol)


def _drgep(net, conds, targets, tols):

    #-----------------------------------
    # Calculate interaction coefficients
    #-----------------------------------

    R_TB = np.zeros(len(net.unique_nuclei), dtype=np.float64)
    adj_nuc = get_adj_nuc(net)

    for comp, rho, T in conds:
        rvals = net.evaluate_rates(rho=rho, T=T, composition=comp)
        _drgep_kernel(net, R_TB, rvals, targets, tols, adj_nuc)

    return R_TB


def _drgep_mpi(net, conds, targets, tols):

    #----------
    # Init. MPI
    #----------

    comm = MPI.COMM_WORLD
    MPI_N = comm.Get_size()
    MPI_rank = comm.Get_rank()

    #-----------------------------------------------
    # Calculate interaction coefficients in parallel
    #-----------------------------------------------

    R_TB_loc = np.zeros(len(net.unique_nuclei), dtype=np.float64)
    adj_nuc = get_adj_nuc(net)

    for i in range(MPI_rank, len(conds), MPI_N):
        comp, rho, T = conds[i]
        rvals = net.evaluate_rates(rho=rho, T=T, composition=comp)
        _drgep_kernel(net, R_TB_loc, rvals, targets, tols, adj_nuc)

    R_TB = np.zeros_like(R_TB_loc)
    comm.Allreduce([R_TB_loc, MPI.DOUBLE], [R_TB, MPI.DOUBLE], op=MPI.MAX)

    return R_TB


def _drgep_numpy(net, conds, targets, tols):

    #----------------------------------------
    # Unpack conditions; make precalculations
    #----------------------------------------

    comp_L, rho_L, T_L = conds

    net.calc_count_matrices()
    net.update_rate_coef_arr()
    adj_nuc = get_adj_nuc(net)

    #------------------------------------------------
    # Calculate interaction coefficients (vectorized)
    #------------------------------------------------

    R_TB = np.zeros(len(net.unique_nuclei), dtype=np.float64)

    for comp in comp_L:
        net.update_yfac_arr(comp)
        for rho in rho_L:
            net.update_prefac_arr(rho, comp)
            for T in T_L:
                rvals_arr = net.evaluate_rates_arr(T)
                _drgep_kernel_numpy(net, R_TB, rvals_arr, targets, tols, adj_nuc)

    net.clear_arrays()
    return R_TB


def _drgep_mpi_numpy(net, conds, targets, tols):

    #----------------------------------------
    # Unpack conditions; make precalculations
    #----------------------------------------

    n = tuple(map(len, conds))
    comp_L, rho_L, T_L = conds

    net.calc_count_matrices()
    net.update_rate_coef_arr()
    adj_nuc = get_adj_nuc(net)

    #-----------------------------------
    # Init. MPI and divide up conditions
    #-----------------------------------

    comm = MPI.COMM_WORLD
    MPI_N = comm.Get_size()
    MPI_rank = comm.Get_rank()

    comp_idx, comp_step, rho_idx, rho_step, T_idx, T_step = mpi_numpy_decomp(MPI_N, MPI_rank, n)

    #--------------------------------------------------------------
    # Calculate interaction coefficients (vectorized and using MPI)
    #--------------------------------------------------------------

    R_TB_loc = np.zeros(len(net.unique_nuclei), dtype=np.float64)

    for i in range(comp_idx, n[0], comp_step):
        comp = comp_L[i]
        net.update_yfac_arr(comp)
        for j in range(rho_idx, n[1], rho_step):
            rho = rho_L[j]
            net.update_prefac_arr(rho, comp)
            for k in range(T_idx, n[2], T_step):
                T = T_L[k]
                rvals_arr = net.evaluate_rates_arr(T)
                _drgep_kernel_numpy(net, R_TB_loc, rvals_arr, targets, tols, adj_nuc)

    R_TB = np.zeros_like(R_TB_loc)
    comm.Allreduce([R_TB_loc, MPI.DOUBLE], [R_TB, MPI.DOUBLE], op=MPI.MAX)

    net.clear_arrays()
    return R_TB


def drgep(net, conds, targets, tols, returnobj='net', use_mpi=False, use_numpy=False):
    """
    Implementation of Directed Relation Graph with Error Propagation (DRGEP) reduction
    method described in Pepiot-Desjardins and Pitch 2008 (doi:10.1016/j.combustflame.2007.10.020)
    and Niemeyer and Sung 2011 (doi:10.1016/j.combustflame.2010.12.010).

    :param net: The network (RateCollection) to reduce.
    :param conds: A set of conditions to reduce over. Should either be a sequence of (composition,
        density, temperature) sequences/tuples if running in standard mode, or a sequence of 3
        sequences ((composition, density, temperature) ordering) if running in NumPy mode. In the
        latter case, the sequences will be permuted to create the dataset. The compositions should
        be pynucastro Composition objects.
    :param targets: A collection of target nuclei (or a single target nucleus) to run the
        graph search algorithm from. Should be supplied as pynucastro Nucleus objects.
    :param tols: Tolerance(s) or cutoff threshold(s) to use for paths from each of the target nuclei.
        Nuclei whose interaction coefficients do not meet the specified tolerance will have their
        interaction coefficients set to 0.0. Can be a single number (will be the same for all targets)
        or a separate value for each target nucleus.
    :param returnobj: The type of object to return. Options are 'net' (a reduced network, the default
        setting), 'nuclei' (unique nuclei with nonzero interaction coefficients, ordered so the
        interaction coefficients are descending), and 'coeff' (the interaction coefficients as a
        NumPy array, with entries corresponding to nuclei in *net.unique_nuclei*).
    :param use_mpi: Whether to divide up the set of conditions across MPI processes or not. Default
        setting is *False*.
    :param use_numpy: Whether to use NumPy to vectorize the interaction coefficient calculations or
        not. This is more memory intensive and may actually hinder performance for some setups.
        Conditions should be supplied as 3 lists that will be permuted to form the dataset (see
        *conds* parameter). Default setting is *False*.
    """

    #------------------
    # Process arguments
    #------------------

    if returnobj not in {'net', 'nuclei', 'coeff'}:
        raise ValueError(f"Invalid 'returnobj' argument: '{returnobj}'.")
    targets = to_list(targets)
    tols = to_list(tols, len(targets))

    #-------------------------------------------------------
    # Determine operation mode and launch appropriate helper
    #-------------------------------------------------------

    if use_mpi:
        if use_numpy:
            R_TB = _drgep_mpi_numpy(net, conds, targets, tols)
        else:
            R_TB = _drgep_mpi(net, conds, targets, tols)
    elif use_numpy:
        R_TB = _drgep_numpy(net, conds, targets, tols)
    else:
        R_TB = _drgep(net, conds, targets, tols)

    #------------------------
    # Return requested object
    #------------------------

    if returnobj == 'net':
        nuclei = [net.unique_nuclei[i] for i in range(len(net.unique_nuclei))
                  if R_TB[i] > 0.0]
        return net.linking_nuclei(nuclei)
    if returnobj == 'nuclei':
        idx = sorted(range(len(R_TB)), key=lambda i: R_TB[i], reverse=True)
        return [net.unique_nuclei[i] for i in idx if R_TB[i] > 0.0]
    if returnobj == 'coeff':
        return R_TB
    # returnobj was checked at the start of the function
    assert False
