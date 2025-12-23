"""Methods to ease the creation of networks."""

from pynucastro.rates import DerivedRate, ReacLibLibrary, TabularLibrary

from .amrexastro_cxx_network import AmrexAstroCxxNetwork
from .fortran_network import FortranNetwork
from .python_network import PythonNetwork
from .simple_cxx_network import SimpleCxxNetwork


def network_helper(nuclei, *,
                   network_type="python",
                   inert_nuclei=None,
                   use_detailed_balance=True,
                   use_tabular_rates=True,
                   tabular_ordering=None,
                   with_reverse=True,
                   verbose=False):
    """Generate a basic network connecting all of the input nuclei
    using all of the rates known to pynucastro.

    Parameters
    ----------
    nuclei : Iterable(Nucleus)
        the nuclei to use for the network
    network_type : str
        The type of network to create.  Allowed values are:

        * "python" : create a :py:obj:`PythonNetwork <pynucastro.networks.python_network.PythonNetwork>`
        * "cxx" : create a :py:obj:`SimpleCxxNetwork <pynucastro.networks.simple_cxx_network.SimpleCxxNetwork>`
        * "fortran" : create a :py:obj:`FortranNetwork <pynucastro.networks.fortran_network.FortranNetwork>`
        * "amrex" : create a :py:obj:`AmrexAstroCxxNetwork <pynucastro.networks.amrexastro_cxx_network.AmrexAstroCxxNetwork>`
    inert_nuclei : list, tuple
        an iterable of Nuclei that should be part of the collection but
        are not linked via reactions to the other Nuclei in the network.
    use_detailed_balanace : bool
        Do we rederive inverse rates using detailed balance?
    use_tabular_rates : bool
        Do we include tabulated weak rates?
    tabular_ordering : Iterable(str)
        If we are including tabular rates, a list of sources
        can be provided to specify which rate sources are used,
        and the priority that each source should have.
    with_reverse : bool
        Do we include the reverse rates from ReacLib?
    verbose : bool
        Output more information

    Returns
    -------
    PythonNetwork, SimpleCxxNetwork, AmrexAstroCxxNetwork, FortranNetwork

    """

    rl = ReacLibLibrary()
    lib = rl.linking_nuclei(nuclei, with_reverse=with_reverse)

    if use_tabular_rates:
        if tabular_ordering:
            tl = TabularLibrary(ordering=tabular_ordering)
        else:
            tl = TabularLibrary()

        lib += tl.linking_nuclei(nuclei)

        # if we have both a tabular and ReacLib rate,
        # remove the ReacLib version

        lib.eliminate_duplicates(rate_type_preference="tabular")

    if use_detailed_balance:
        rates_to_derive = lib.backward().get_rates()

        # now for each of those derived rates, look to see if the pair exists
        for r in rates_to_derive:
            fr = lib.get_rate_by_nuclei(r.products, r.reactants)
            if fr:
                if verbose:
                    print(f"modifying {r} from {fr}")
                lib.remove_rate(r)
                d = DerivedRate(source_rate=fr, use_pf=True, use_unreliable_spins=True)
                lib.add_rate(d)

    if network_type == "python":
        return PythonNetwork(libraries=[lib],
                             inert_nuclei=inert_nuclei, verbose=verbose)

    if network_type == "cxx":
        return SimpleCxxNetwork(libraries=[lib],
                                inert_nuclei=inert_nuclei, verbose=verbose)

    if network_type == "fortran":
        return FortranNetwork(libraries=[lib],
                              inert_nuclei=inert_nuclei, verbose=verbose)

    if network_type == "amrex":
        return AmrexAstroCxxNetwork(libraries=[lib],
                                    inert_nuclei=inert_nuclei, verbose=verbose)

    raise ValueError("invalid network_type")
