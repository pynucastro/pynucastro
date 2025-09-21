#!/usr/bin/env python3

import argparse
import sys
import time
from collections import namedtuple

import numpy as np

from pynucastro import Nucleus
from pynucastro.reduction import mpi_importer, sens_analysis
from pynucastro.reduction.drgep_module import drgep
from generate_data import dataset
from load_network import load_network

MPI = mpi_importer()

NetInfo = namedtuple("NetInfo", "y ydot z a ebind m")


def rel_err(x, x0):
    """Compute the relative error between two NumPy arrays."""
    return np.abs((x - x0) / x0)


def get_errfunc_enuc(net_old, conds):
    """Compute the error in nuclear energy generation."""

    enucdot_list = []

    for comp, rho, T in conds:
        enucdot_list.append(net_old.evaluate_energy_generation(rho, T, comp))

    def erf(net_new):

        err = 0.0
        for (comp, rho, T), enucdot_old in zip(conds, enucdot_list):
            enucdot_new = net_new.evaluate_energy_generation(rho, T, comp)
            err = max(err, rel_err(enucdot_new, enucdot_old))

        return err

    return erf


def main():
    """Initialize and run the reduction."""

    #-----------------------------------
    # Setup parser and process arguments
    #-----------------------------------

    description = "Example/test script for running a selected reduction algorithm."
    algorithm_help = "The algorithm to use. Currently only supports 'drgep'."
    endpoint_help = "The nucleus to use as an endpoint in the unreduced network."
    datadim_help = """The dimensions of the dataset to generate. Ordering is the number of densities,
            then temperatures, then metallicities. Default is [4, 4, 4]."""
    drho_help = "Number of density points to include in the dataset."
    dtemp_help = "Number of temperature points to include in the dataset."
    dmetal_help = "Number of metallicity points to include in the dataset."
    brho_help = "Range of densities to include in the dataset (in g/cm^3)."
    btemp_help = "Range of temperatures to include in the dataset (in K)."
    bmetal_help = """Range of metallicities to include in the dataset. Half of the metallicity will
            be in C12, and the rest will be divided evenly among the remaining nuclei."""
    library_help = "Name of the library to use to load the network."
    targets_help = """Target nuclei to use for the algorithm. Will be protons and the endpoint by
            default."""
    tol_help = """Tolerance(s) to use for the algorithm. Will be [1e-3] + [1e-2]*(len(targets)-1)
            by default for 'drgep'."""
    use_mpi_help = "Enable MPI for this run."
    use_numpy_help = """If the algorithm has a 'use_numpy' option, turn it on. Some algorithms always
            run as if use_numpy=True."""
    sa_help = "Error threshold to use when performing sensitivity analysis."

    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-a', '--algorithm', default='drgep', help=algorithm_help)
    parser.add_argument('-e', '--endpoint', default=Nucleus('te108'), type=Nucleus,
            help=endpoint_help)
    parser.add_argument('-d', '--datadim', nargs=3, default=[4]*3, help=datadim_help)
    parser.add_argument('-dr', '--drho', type=int, help=drho_help)
    parser.add_argument('-dt', '--dtemp', help=dtemp_help)
    parser.add_argument('-dz', '--dmetal', type=int, help=dmetal_help)
    parser.add_argument('-br', '--brho', nargs=2, type=float, help=brho_help)
    parser.add_argument('-bt', '--btemp', nargs=2, type=float, help=btemp_help)
    parser.add_argument('-bz', '--bmetal', nargs=2, type=float, help=bmetal_help)
    parser.add_argument('-l', '--library', default="rp-process-lib", help=library_help)
    parser.add_argument('-t', '--targets', nargs='*', type=Nucleus, help=targets_help)
    parser.add_argument('--tol', nargs='*', type=float, help=tol_help)
    parser.add_argument('--use_mpi', action='store_true', help=use_mpi_help)
    parser.add_argument('--use_numpy', action='store_true', help=use_numpy_help)
    parser.add_argument('-s', '--sens_analysis', default=0.05, type=float, help=sa_help)
    args = parser.parse_args(sys.argv[1:])

    args.algorithm = args.algorithm.lower()

    if args.algorithm == 'drgep':
        alg = drgep
        permute = not args.use_numpy
    else:
        raise ValueError(f'Algorithm {args.algorithm} is not supported.')

    if args.drho:
        args.datadim[0] = args.drho
    if args.dtemp:
        args.datadim[1] = args.dtemp
    if args.dmetal:
        args.datadim[2] = args.dmetal

    if not args.targets:
        args.targets = [Nucleus('p'), args.endpoint]

    if not args.tol:
        if args.algorithm == 'drgep':
            args.tol = [1e-3] + [1e-2] * (len(args.targets)-1)
    else:
        if args.algorithm == 'drgep':
            if len(args.tol) == 1:
                args.tol = [args.tol] * len(args.targets)
            elif len(args.tol) != len(args.targets):
                raise ValueError(
                    f"For '{args.algorithm}', there should be one tolerance for"
                    " each target nucleus."
                )

    #-----------------------------
    # Load the network and dataset
    #-----------------------------

    net = load_network(args.endpoint, library_name=args.library)
    data = list(dataset(net, args.datadim, permute, args.brho, args.btemp, args.bmetal))

    #-----------------------------
    # Prepare to run the algorithm
    #-----------------------------

    alg_args = \
    {
        'net': net,
        'conds': data,
        'targets': args.targets,
        'tols': args.tol,
        'returnobj': 'nuclei',
        'use_mpi': args.use_mpi,
        'use_numpy': args.use_numpy
    }

    if args.algorithm == 'drgep':
        args.algorithm = args.algorithm.upper()

    #------------------
    # Run the algorithm
    #------------------

    if args.use_mpi and MPI.COMM_WORLD.Get_rank() == 0:
        print(f"Commencing reduction with {MPI.COMM_WORLD.Get_size()} processes.\n")
    elif not args.use_mpi:
        print("Commencing reduction without MPI.\n")

    t0 = time.time()
    nuclei = alg(**alg_args)
    dt = time.time() - t0

    red_net = net.linking_nuclei(nuclei)
    second_data = list(dataset(net, args.datadim, True, args.brho, args.btemp, args.bmetal))
    errfunc = get_errfunc_enuc(net, second_data)

    if not args.use_mpi or MPI.COMM_WORLD.Get_rank() == 0:

        print()
        print(f"{args.algorithm} reduction took {dt:.3f} s.")
        print("Number of species in full network: ", len(net.unique_nuclei))
        print(f"Number of species in {args.algorithm} reduced network: ", len(nuclei))
        print("Reduced Network Error:", f"{errfunc(red_net)*100:.2f}%")
        print()

    # Perform sensitivity analysis
    if args.use_mpi:
        MPI.COMM_WORLD.Barrier()
    t0 = time.time()
    red_net = sens_analysis(red_net, errfunc, args.sens_analysis, args.use_mpi)
    dt = time.time() - t0

    if not args.use_mpi or MPI.COMM_WORLD.Get_rank() == 0:

        print(f"Greedy sensitivity analysis reduction took {dt:.3f} s.")
        print(f"Number of species in {args.algorithm} + sensitivity analysis reduced network: ",
                len(red_net.unique_nuclei))
        print("Error: ", f"{errfunc(red_net)*100:.2f}%")

    print("final network:")
    print(red_net.summary())


if __name__ == "__main__":
    main()
