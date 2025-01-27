# XRB rate module generator
from pynucastro import ReacLibLibrary, PythonNetwork
from pynucastro.constants import constants
from pynucastro.nucdata import Nucleus
from pynucastro.reduction import mpi_importer, sens_analysis
from pynucastro.reduction.drgep import drgep
from pynucastro.reduction.generate_data import dataset

from rp_process import rp_process
from reduction import get_net_info, enuc_dot, rel_err, get_errfunc_enuc

def make_rp_process(endpoint):
    ############################################
    # Create rp-process network up to endpoint #
    ############################################
    return rp_process(endpoint)

def make_reduced(net, use_MPI, use_numpy, alg_tol=None, sens_tol = 0.05, 
        MPI=None):
    ##########################
    # Reduce created network #
    ##########################
    # Make sure MPI makes sense
    assert use_MPI and MPI
    # Setup the conditions
    Nrho = 4
    Ntemp = 4
    Nmetal = 4
    datadim = [Nrho, Ntemp, Nmetal]
    permute = not use_numpy #Lets algorithm use numpy
    data1 = list(dataset(net, datadim, permute))
    data2 = list(dataset(net, datadim, True))
    # Determine targets
    targets = [Nucleus('p'), net.get_nuclei()[-1]]
    # Create the error function and set tolerences
    errfunc = get_errfunc_enuc(net, data2)
    if not alg_tol:
        alg_tol = [1e-3] + [1e-2] * (len(targets)-1)
    else:
        if type(alg_tol) is float:
            alg_tol = [alg_tol] * len(targets)
        if len(alg_tol) == 1:
            alg_tol = [alg_tol] * len(targets)
        elif len(alg_tol) != len(targets):
            raise ValueError(
                    f"For 'drgep', there hsould be one tolerance for each"
                    " target nucleus"
                    )
    # Construct the algorithm arguments
    alg_args = \
        {
                'net': net,
                'conds': data1,
                'targets': targets,
                'tols': alg_tol,
                'returnobj': 'nuclei',
                'use_mpi': use_MPI,
                'use_numpy': use_numpy
        }
    # Run the algorithm
    red_nuclei = drgep(**alg_args)
    reduced_net = net.linking_nuclei(red_nuclei)
    # Perform sensitivity analysis
    if use_MPI:
        MPI.COMM_WORLD.Barrier()
    reduced_net = sens_analysis(reduced_net, errfunc, sens_tol, use_MPI)

    return reduced_net

def make_CNO():
    ############################
    # Create CNO extra network #
    ############################
    rl = ReacLibLibrary()
    h_burn = rl.linking_nuclei(["h1", "he4",
                                "c12", "c13",
                                "n13", "n14", "n15",
                                "o14", "o15", "o16","o17","o18",
                                "f17", "f18","f19",
                                "ne18", "ne19", "ne20",
                                "mg22", "mg24"],
                               with_reverse=False)

    return PythonNetwork(libraries=[h_burn], inert_nuclei=["fe56"])

def main():
    # Check if MPI can be used and use if it can
    MPI = mpi_importer()
    use_MPI = (MPI.COMM_WORLD.Get_size() > 1)

    endpoint = Nucleus('ni56')

    rp_net = make_rp_process(endpoint)
    reduced_net = make_reduced(rp_net, use_MPI, False, None, 0.05, MPI)
    CNO_net = make_CNO()

    # The rest should only be done by one process
    if not use_MPI or MPI.COMM_WORLD.Get_rank() == 0:
        # See difference between reduced and original net
        print(f"The original network has {len(rp_net.get_nuclei())} nuclei")
        print(f"The reduced network has {len(reduced_net.get_nuclei())} nuclei")

        # Write the networks
        rp_net.write_network("rp"+str(endpoint)+"_rhs.py")
        reduced_net.write_network("reduced_rhs.py")
        CNO_net.write_network("CNO_rhs.py")

        # Print the networks
        nets = [rp_net, reduced_net, CNO_net]
        fnames = ["rp"+str(endpoint), "reduced", "CNOextra"]

        for net, fname in zip(nets, fnames):
            net.plot(outfile = fname + ".png")

if __name__ == "__main__":
    main()
