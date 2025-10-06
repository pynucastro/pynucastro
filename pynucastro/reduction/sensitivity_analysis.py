"""Methods and functions used for sensitivity analysis."""

import numpy as np

from pynucastro.mpi_utils import mpi_importer
from pynucastro.nucdata import Nucleus

MPI = mpi_importer()


def binary_search_trim(network, nuclei, errfunc, *, thresh=0.05, args=None):
    """Given an array of nuclei sorted roughly by relative importance,
    perform a binary search to trim out nuclei from the network until
    the error is as close as possible to the given threshold without
    exceeding it. Nuclei whose removal will result in only a small
    increase in error need to be packed towards the back of the array
    for the binary search to work effectively.

    Parameters
    ----------
    network : RateCollection
        The network to reduce.
    nuclei : Iterable(Nucleus) or Iterable(str)
        Nuclei to consider for the final network, sorted by decreasing
        importance (i.e.  most important nuclei first). Importance can
        be determined by something like the *drgep* function.
    errfunc : Callable
        Error function to use when evaluating error, with the signature
        ``error(net, *args)``, where ``net`` is the reduced network as
        an argument and return the relative error produced by the reduction.
        If ``use_mpi`` is ``False``, the error function can be parallelized
        with MPI. Otherwise ``sens_analysis`` will be parallelized and the
        error function should not be.
    thresh : float
        Threshold for acceptable error. Default is 0.05.
    args : tuple
        Additional arguments to pass through to the error function

    Returns
    -------
    net : RateCollection
        A reduced reaction network with an evaluated error
        approximately equal to the supplied threshold.

    """

    nuclei = Nucleus.cast_list(nuclei)

    start_idx = 0
    seg_size = len(nuclei) / 2

    while seg_size >= 0.5:

        # Divide up into segments
        end_idx = start_idx + round(seg_size)
        red_net = network.linking_nuclei(nuclei[:end_idx])

        # Evaluate error
        if args is None:
            err = errfunc(red_net)
        else:
            err = errfunc(red_net, *args)

        if err <= thresh:
            seg_size /= 2
        else:
            start_idx += round(seg_size)
            seg_size /= 2

    return network.linking_nuclei(nuclei[:start_idx+1])


def _progress_bar(frac, size=50):

    n = round(size*frac)
    progress_bar = '[' + '⊙'*n + ' '*(size-n) + ']'
    if frac < 1.0:
        end = '\r'
    else:
        end = '\n'
    print(progress_bar, f'{round(100*frac)}%', end=end)


def sens_analysis(network, errfunc, *, thresh=0.05, args=None,
                  use_mpi=False, print_prog=True):
    """Given a reaction network, remove nuclei from the network
    one-by-one until the induced error is as close to the given
    threshold as possible without exceeding it. This will test nuclei
    for removal individually and remove the one that induces the
    smallest error on each pass. Since it requires O(n^2) error
    function evaluations, this routine is much more expensive than
    ``binary_search``, but it will typically trim the network down
    significantly more.

    Parameters
    ----------
    network : RateCollection
        The network to reduce.  Can be a RateCollection or a subclass.
    errfunc : Callable
        Error function to use when evaluating error, with the signature
        ``error(net, *args)``, where ``net`` is the reduced network as
        an argument and return the relative error produced by the reduction.
        If ``use_mpi`` is ``False``, the error function can be parallelized
        with MPI. Otherwise ``sens_analysis`` will be parallelized and the
        error function should not be.
    thresh : float
        Threshold for acceptable error. Default is 0.05.
    args : tuple
        Additional arguments to pass through to the error function
    use_mpi : bool
        Whether to parallelize the loop over nuclei with each pass or
        not using MPI.  For p MPI processes, the parallelized
        function will require O(n^2/p) error function evaluations per
        process. This option is ``False`` by default. If the error
        function is parallelized using MPI, this option should be set
        to ``False``.
    print_prog : bool
        Whether to print out the progress of the function as it runs
        or not. Includes a progress bar for each pass and messages
        indicating when the algorithm starts and ends.

    Returns
    -------
    net : RateCollection
        A reduced reaction network with an evaluated error
        approximately equal to the supplied threshold.

    """

    if use_mpi:
        comm = MPI.COMM_WORLD
        MPI_N = comm.Get_size()
        MPI_rank = comm.Get_rank()
    else:
        MPI_N = 1
        MPI_rank = 0

    nuclei = list(network.unique_nuclei)
    err = 0.0

    nrem = 0
    print_prog = print_prog and (MPI_rank == 0)

    if print_prog:
        print("Performing sensitivity analysis...")

    while True:

        err = float('inf')

        for i in range(MPI_rank, len(nuclei), MPI_N):

            if print_prog:
                print(f"Pass {nrem+1}:", end=' ')
                _progress_bar(i/len(nuclei))

            nuc = nuclei.pop(i)
            if args is None:
                err_i = errfunc(network.linking_nuclei(nuclei, print_warning=False))
            else:
                err_i = errfunc(network.linking_nuclei(nuclei, print_warning=False), *args)

            if err_i < err:
                err = err_i
                min_idx = i
            nuclei.insert(i, nuc)

        if use_mpi:

            err = comm.gather(err, root=0)
            min_idx = comm.gather(min_idx, root=0)

            if MPI_rank == 0:
                min_rank = np.argmin(err)
                err = err[min_rank]
                min_idx = min_idx[min_rank]

            err = comm.bcast(err, root=0)
            min_idx = comm.bcast(min_idx, root=0)

        if print_prog:
            print(f"Pass {nrem+1}:", end=' ')
            _progress_bar(1.0)

        if err <= thresh:
            nuclei.pop(min_idx)
            nrem += 1
        else:
            break

    if print_prog:
        print(f"Done. Removed {nrem} nuclei.")
    return network.linking_nuclei(nuclei)
